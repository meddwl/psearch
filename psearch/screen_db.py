#!/usr/bin/env python3
# author          : Pavel Polishchuk
# date            : 23.08.2019
# license         : BSD-3
# ==============================================================================

import os
import argparse
from collections import namedtuple
from pmapper.pharmacophore import Pharmacophore
from multiprocessing import Pool
from functools import partial
from rdkit import Chem
from rdkit.Chem import AllChem
from psearch.database import DB


path_query = os.path.join(os.path.split(os.path.abspath(os.path.dirname(__file__)))[0], 'pharmacophores', 'pharms_chembl')
Model = namedtuple('Model', ['name', 'fp', 'pharmacophore', 'output_filename'])
Conformer = namedtuple('Conformer', ['stereo_id', 'conf_id', 'fp', 'pharmacophore'])


def create_parser():
    parser = argparse.ArgumentParser(description='Screen DB with compounds against pharmacophore queries.',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-d', '--dbname', metavar='FILENAME.dat', type=str, required=True,
                        help='input database with generated conformers and pharmacophores.')
    parser.add_argument('-q', '--query', metavar='FILENAME(S) or DIRNAME(S)', type=str, nargs='+',
                        default=[os.path.join(path_query, q) for q in os.listdir(path_query)],
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
                        help='specify if sdf output with matched conformers is required. These files will be created '
                             'in the place as text files.')
    parser.add_argument('--conf', action='store_true', default=False,
                        help='return all conformers matches as separate hits in a hit list. Required to calculate the '
                             'score by the conformer coverage approach (CCA).')
    parser.add_argument('-c', '--ncpu', metavar='INTEGER', default=1, type=int,
                        help='number of cores to use. Default: 1.')
    return parser


def load_confs(mol_name, db):
    bin_step = db.get_bin_step()
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


def screen(mol_name, db, models, output_sdf, match_first_conf):

    def compare_fp(query_fp, fp):
        return (query_fp & fp) == query_fp

    get_transform_matrix = output_sdf

    confs = load_confs(mol_name, db)

    output = []
    for model in models:
        for conf in confs:
            if compare_fp(model.fp, conf.fp):
                res = conf.pharmacophore.fit_model(model.pharmacophore, get_transform_matrix=get_transform_matrix)
                if res:
                    if get_transform_matrix:
                        output.append((mol_name, conf.stereo_id, conf.conf_id, model.output_filename, res[1]))
                    else:
                        output.append((mol_name, conf.stereo_id, conf.conf_id, model.output_filename))
                    if match_first_conf:
                        break
    return output


def save_results(results, output_sdf, db):
    for items in results:
        mol_name, stereo_id, conf_id, out_fname = items[:4]
        if not os.path.exists(os.path.dirname(out_fname)):
            os.makedirs(os.path.dirname(out_fname))
        with open(out_fname, 'at') as f:
            f.write('\t'.join((mol_name, str(stereo_id), str(conf_id))) + '\n')
    if output_sdf:
        for mol_name, stereo_id, conf_id, out_fname, matrix in results:
            m = db.get_mol(mol_name)[stereo_id]
            AllChem.TransformMol(m, matrix, conf_id)
            m.SetProp('_Name', f'{mol_name}-{stereo_id}-{conf_id}')
            with open(os.path.splitext(out_fname)[0] + '.sdf', 'a') as f:
                w = Chem.SDWriter(f)
                w.write(m)
                w.close()


def screen_db(db_fname, queries, output, output_sdf, match_first_conf, min_features, ncpu):

    if output.endswith('.txt') or output.endswith('.sdf'):
        if not os.path.exists(os.path.dirname(output)):
            os.makedirs(os.path.dirname(output), exist_ok=True)
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

    # print(comp_names)

    if ncpu == 1:
        for comp_name in comp_names:
            res = screen(mol_name=comp_name, db=db, models=models, output_sdf=output_sdf, match_first_conf=match_first_conf)
            if res:
                save_results(res, output_sdf, db)
    else:
        p = Pool(ncpu)
        for res in p.imap_unordered(partial(screen, db=db, models=models, output_sdf=output_sdf, match_first_conf=match_first_conf), comp_names, chunksize=10):
            if res:
                save_results(res, output_sdf, db)
        p.close()

    # remove output dir if it is empty
    # if os.path.exists(output) and os.path.isdir(output) and not os.listdir(output):
    #     os.rmdir(output)


def entry_point():
    parser = create_parser()
    args = parser.parse_args()
    screen_db(db_fname=args.dbname,
              queries=args.query,
              output=args.output,
              output_sdf=args.output_sdf,
              match_first_conf=not args.conf,
              min_features=args.min_features,
              ncpu=args.ncpu)


if __name__ == '__main__':
    entry_point()
