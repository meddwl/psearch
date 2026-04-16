#!/usr/bin/env python3
# author          : Pavel Polishchuk
# date            : 23.08.2019
# license         : BSD-3
# ==============================================================================

import os
import sys
import time
import argparse
from collections import namedtuple
from pmapper.pharmacophore import Pharmacophore
from multiprocessing import Pool
from functools import partial
from rdkit import Chem
from rdkit.Chem import AllChem
from psearch.database import DB


path_query = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'pharmacophores', 'chembl_models')
Model = namedtuple('Model', ['name', 'fp', 'pharmacophore', 'output_filename'])
Conformer = namedtuple('Conformer', ['stereo_id', 'conf_id', 'fp', 'pharmacophore'])


def create_parser():
    """Build the CLI argument parser for screen_db."""
    parser = argparse.ArgumentParser(description='Screen DB with compounds against pharmacophore queries.',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-d', '--dbname', metavar='FILENAME.dat', type=str, required=True,
                        help='input database with generated conformers and pharmacophores.')
    parser.add_argument('-q', '--query', metavar='FILENAME(S) or DIRNAME(S)', type=str, nargs='+', default=None,
                        help='pharmacophore model(s) or directory path(s). If a directory is specified all '
                             'pma- and xyz-files will be used for screening as pharmacophore models.'
                             'The ligand-based pharmacophore models, that created from the ChEMBL database '
                             'using  the psearch tool, are used by default.')
    parser.add_argument('-o', '--output', metavar='FILENAME or DIRNAME', required=True, type=str,
                        help='a text (.txt) file which will store names of compounds which fit the model. In the case '
                             'multiple query models or directories were supplied as input'
                             'this should be the path to a directory where output files will be created to store '
                             'screening results. If multiple directories were specified as input the corresponding '
                             'directories will be created in the output dir. Names of created  directories will be '
                             'taken from the bottom level of input directories, e.g. path/to/model/ will be stored in '
                             'output_dir/model. Beware, existed output files/directories will be overwritten.')
    parser.add_argument('-f', '--min_features', metavar='INTEGER', default=None, type=int,
                        help='minimum number of features with distinct coordinates in models. Models having less '
                             'number of features will be skipped. Default: all models will be screened.')
    parser.add_argument('-z', '--output_sdf', action='store_true', default=False,
                        help='Write an SDF file alongside each hit-list file containing the matching 3D conformers.')
    parser.add_argument('--conf', action='store_true', default=False,
                        help='Report each matching conformer as a separate hit. '
                             'Required for conformer-coverage-approach (CCA) scoring.')
    parser.add_argument('-c', '--ncpu', metavar='INTEGER', default=1, type=int,
                        help='number of cores to use. Default: 1.')
    parser.add_argument('-v', '--verbose', action='store_true', default=False,
                        help='print progress to STDERR.')
    return parser


def load_confs(mol_name, db, bin_step):
    """Load all conformers for a compound from the database and reconstruct Pharmacophore objects.

    Args:
        mol_name: Compound identifier string.
        db: Open DB instance.
        bin_step: Bin width (Å) used to build the Pharmacophore objects; must match the database.

    Returns:
        List of Conformer namedtuples (stereo_id, conf_id, fp, pharmacophore).
    """
    fp_dict = db.get_fp(mol_name)
    ph_dict = db.get_pharm(mol_name)
    res = []
    for stereo_id in fp_dict:
        try:
            for conf_id, (fp, coord) in enumerate(zip(fp_dict[stereo_id], ph_dict[stereo_id])):
                p = Pharmacophore(bin_step=bin_step)
                p.load_from_feature_coords(coord)
                res.append(Conformer(stereo_id, conf_id, fp, p))
        except:
            print(mol_name)
    return res


def read_models(queries, output, bin_step, min_features):
    """Parse pharmacophore model files and prepare them for screening.

    Accepts either a list of .pma/.xyz files or a list of directories. Models with
    fewer than `min_features` distinct-coordinate features are skipped.

    Args:
        queries: List of file paths or directory paths containing .pma/.xyz model files.
        output: Output path (file or directory) used to derive per-model output filenames.
        bin_step: Bin width (Å) used to update pharmacophore fingerprints.
        min_features: Skip models with fewer distinct-coordinate features than this value.
                      None screens all models.

    Returns:
        List of Model namedtuples (name, fp, pharmacophore, output_filename).
    """
    if all(os.path.isdir(item) for item in queries):
        input_fnames = []
        output_fnames = []
        for dname in queries:
            dname = os.path.abspath(dname)
            for f in os.listdir(dname):
                if os.path.isfile(os.path.join(dname, f)) and (f.endswith('.pma') or f.endswith('.xyz')):
                    input_fnames.append(os.path.join(dname, f))
                    output_fnames.append(os.path.join(output, os.path.basename(dname), os.path.splitext(os.path.basename(f))[0] + '.txt'))
    elif all(os.path.isfile(item) for item in queries):
        input_fnames = tuple(os.path.abspath(f) for f in queries if f.endswith('.pma') or f.endswith('.xyz'))
        output_fnames = tuple(os.path.join(output, os.path.splitext(os.path.basename(f))[0] + '.txt') for f in input_fnames)
    else:
        raise ValueError('Input queries should be all either files or directories not a mix.')

    res = []
    for input_fname, output_fname in zip(input_fnames, output_fnames):
        p = Pharmacophore()
        if input_fname.endswith('.pma'):
            p.load_from_pma(input_fname)
        elif input_fname.endswith('.xyz'):
            p.load_from_xyz(input_fname)
        # skip models with less number of features with distinct coordinates that given
        if min_features is not None and len(set(xyz for label, xyz in p.get_feature_coords())) < min_features:
            continue
        p.update(bin_step=bin_step)
        fp = p.get_fp()
        res.append(Model(input_fname, fp, p, output_fname))

    return res


def screen(mol_name, db, models, output_sdf, match_first_conf, bin_step):
    """Screen all conformers of a compound against a list of pharmacophore models.

    For each model, iterates over conformers and uses a fingerprint pre-filter before
    full pharmacophore fitting. Returns only matching (model, conformer) combinations.

    Args:
        mol_name: Compound identifier string.
        db: Open DB instance.
        models: List of Model namedtuples from read_models.
        output_sdf: If True, also retrieve the transformation matrix and RMSD for SDF output.
        match_first_conf: If True, stop after the first matching conformer per model
                          (faster; use False for CCA scoring).
        bin_step: Bin width (Å) used to load conformers; must match the database.

    Returns:
        List of tuples (mol_name, stereo_id, conf_id, output_filename) or
        (mol_name, stereo_id, conf_id, output_filename, matrix, rms) when output_sdf is True.
    """
    def compare_fp(query_fp, fp):
        """Return True if query_fp is a subset of fp (all query bits are set in the molecule fp)."""
        return (query_fp & fp) == query_fp

    get_transform_matrix = output_sdf
    get_rms = output_sdf

    confs = load_confs(mol_name, db, bin_step)

    output = []
    for model in models:
        for conf in confs:
            if compare_fp(model.fp, conf.fp):
                res = conf.pharmacophore.fit_model(model.pharmacophore,
                                                   get_transform_matrix=get_transform_matrix, get_rms=get_rms)
                if res:
                    if get_transform_matrix:
                        output.append((mol_name, conf.stereo_id, conf.conf_id, model.output_filename, res[1], res[2]))
                    else:
                        output.append((mol_name, conf.stereo_id, conf.conf_id, model.output_filename))
                    if match_first_conf:
                        break
    return output


def save_results(results, output_sdf, db):
    """Write screening hits to text hit-list files and, optionally, to SDF files.

    Hit-list files contain one tab-separated line per hit: mol_name, stereo_id, conf_id.
    SDF files (written when output_sdf is True) contain the matching 3D conformer
    superimposed onto the pharmacophore model, with the RMSD stored as an SD property.

    Args:
        results: List of tuples returned by screen — either 4-element (text only) or
                 6-element (text + matrix + rms) when output_sdf is True.
        output_sdf: If True, also write SDF files alongside the hit-list files.
        db: Open DB instance (needed to retrieve 3D coordinates for SDF output).
    """
    # Group by output filename to batch writes and minimise open/close calls
    created_dirs = set()
    by_fname = {}
    for items in results:
        mol_name, stereo_id, conf_id, out_fname = items[:4]
        by_fname.setdefault(out_fname, []).append(items)

    for out_fname, items_list in by_fname.items():
        out_dir = os.path.dirname(out_fname)
        if out_dir and out_dir not in created_dirs:
            os.makedirs(out_dir, exist_ok=True)
            created_dirs.add(out_dir)
        with open(out_fname, 'at') as f:
            for items in items_list:
                mol_name, stereo_id, conf_id = items[0], items[1], items[2]
                f.write('\t'.join((mol_name, str(stereo_id), str(conf_id))) + '\n')

    if output_sdf:
        sdf_by_fname = {}
        for mol_name, stereo_id, conf_id, out_fname, matrix, rms in results:
            sdf_by_fname.setdefault(out_fname, []).append((mol_name, stereo_id, conf_id, matrix, rms))
        for out_fname, sdf_items in sdf_by_fname.items():
            with open(os.path.splitext(out_fname)[0] + '.sdf', 'a') as f:
                w = Chem.SDWriter(f)
                for mol_name, stereo_id, conf_id, matrix, rms in sdf_items:
                    m = db.get_mol(mol_name)[stereo_id]
                    AllChem.TransformMol(m, matrix, conf_id)
                    m.SetProp('_Name', f'{mol_name}-{stereo_id}-{conf_id}')
                    m.SetProp("RMSD", str(round(rms, 4)))
                    w.write(m)
                w.close()


def screen_db(db_fname, queries, output, output_sdf, match_first_conf, min_features, ncpu, verbose):
    """Orchestrate parallel pharmacophore virtual screening of a psearch database.

    Reads all compound conformers from `db_fname`, screens them against `queries`,
    and writes hit lists (and optionally SDF files) to `output`.

    Args:
        db_fname: Path to the psearch database (.dat file).
        queries: List of .pma/.xyz file paths or directory paths containing model files.
        output: Output file path (single model) or directory path (multiple models).
        output_sdf: If True, write matched conformers to SDF files alongside hit lists.
        match_first_conf: If True, stop after the first conformer match per model (faster).
                          Set to False when all matching conformers are needed (e.g. CCA).
        min_features: Skip models with fewer distinct-coordinate features. None screens all.
        ncpu: Number of parallel worker processes.
        verbose: If True, print progress to stderr every 10 molecules.
    """
    start_time = time.time()

    if output.endswith('.txt') or output.endswith('.sdf'):
        output_dir = os.path.dirname(os.path.abspath(output))
        if not os.path.exists(output_dir):
            os.makedirs(output_dir, exist_ok=True)
    else:
        if not os.path.exists(output):
            os.makedirs(output, exist_ok=True)

    if output.endswith('.sdf'):  # forcibly set output format
        output_sdf = True

    db = DB(db_fname, flag='r')
    bin_step = db.get_bin_step()
    models = read_models(queries, output, bin_step, min_features)   # return list of Model namedtuples
    for model in models:
        if os.path.isfile(model.output_filename):
            os.remove(model.output_filename)
        if output_sdf and os.path.isfile(os.path.splitext(model.output_filename)[0] + '.sdf'):
            os.remove(os.path.splitext(model.output_filename)[0] + '.sdf')

    comp_names = db.get_mol_names()

    if ncpu == 1:
        for i, comp_name in enumerate(comp_names, 1):
            res = screen(mol_name=comp_name, db=db, models=models, output_sdf=output_sdf,
                         match_first_conf=match_first_conf, bin_step=bin_step)
            if res:
                save_results(res, output_sdf, db)
            if verbose and i % 10 == 0:
                current_time = time.strftime("%H:%M:%S", time.gmtime(time.time() - start_time))
                sys.stderr.write('\r{} molecules passed/conformers {}'.format(i, current_time))
                sys.stderr.flush()
    else:
        with Pool(ncpu) as p:
            for i, res in enumerate(p.imap_unordered(partial(screen, db=db, models=models, output_sdf=output_sdf,
                                                match_first_conf=match_first_conf, bin_step=bin_step),
                                                comp_names, chunksize=10), 1):
                if res:
                    save_results(res, output_sdf, db)
                if verbose and i % 10 == 0:
                    current_time = time.strftime("%H:%M:%S", time.gmtime(time.time() - start_time))
                    sys.stderr.write('\r{} molecules screened {}'.format(i, current_time))
                    sys.stderr.flush()

    # remove output dir if it is empty
    # if os.path.exists(output) and os.path.isdir(output) and not os.listdir(output):
    #     os.rmdir(output)


def entry_point():
    """CLI entry point for screen_db: parse arguments and call screen_db."""
    parser = create_parser()
    args = parser.parse_args()
    screen_db(db_fname=args.dbname,
              queries=[os.path.join(path_query, q) for q in os.listdir(path_query)] if not args.query else args.query,
              output=args.output,
              output_sdf=args.output_sdf,
              match_first_conf=not args.conf,
              min_features=args.min_features,
              ncpu=args.ncpu,
              verbose=args.verbose)


if __name__ == '__main__':
    entry_point()
