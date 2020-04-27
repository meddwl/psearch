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

Model = namedtuple('Model', ['name', 'fp', 'pharmacophore', 'output_filename'])
Conformer = namedtuple('Conformer', ['stereo_id', 'conf_id', 'fp', 'pharmacophore'])


def create_parser():
    parser = argparse.ArgumentParser(description='Screen SQLite DB with compounds against pharmacophore queries.',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-d', '--dbname', metavar='FILENAME', type=str, required=True,
                        help='input database with generated conformers and pharmacophores.')
    parser.add_argument('-q', '--query', metavar='FILENAME(S) or DIRNAME', required=True, type=str, nargs='+',
                        help='pharmacophore model or models or a directory path. If a directory is specified all '
                             'pma- and xyz-files will be used for screening as pharmacophore models.')
    parser.add_argument('-o', '--output', metavar='FILENAME or DIRNAME', required=True, type=str,
                        help='a text (.txt) file which will store names of compounds which fit the model or '
                             'a sdf file which will store matched conformers of compounds. The output format '
                             'will be recognized by file extension. In the case multiple query models were supplied '
                             'this should be the path to a directory where output files will be created to store '
                             'screening results. In this case file format should be specified by a separate argument. '
                             'Existed output files will be overwritten.')
    parser.add_argument('-f', '--min_features', metavar='INTEGER', default=None, type=int,
                        help='minimum number of features with distinct coordinates in models. Models having less '
                             'number of features will be skipped. Default: all models will be screened.')
    parser.add_argument('-z', '--output_sdf', action='store_true', default=False,
                        help='specify if sdf output with matched conformers is required.')
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
        for conf_id, (fp, coord) in enumerate(zip(fp_dict[stereo_id], ph_dict[stereo_id])):
            p = Pharmacophore(bin_step=bin_step)
            p.load_from_feature_coords(coord)
            res.append(Conformer(stereo_id, conf_id, fp, p))
    return res


def read_models(queries, output, is_output_sdf, bin_step, min_features):

    if len(queries) == 1 and os.path.isdir(queries[0]):
        input_fnames = tuple(os.path.abspath(os.path.join(queries[0], f)) for f in os.listdir(queries[0]) if f.endswith('.pma') or f.endswith('.xyz'))
    else:
        input_fnames = tuple(os.path.abspath(f) for f in queries if os.path.isfile(f) and f.endswith('.pma'))

    model_names = tuple(os.path.splitext(os.path.basename(f))[0] for f in input_fnames)

    output = os.path.abspath(output)
    if os.path.isdir(output):
        ext = '.sdf' if is_output_sdf else '.txt'
        output_fnames = tuple(os.path.join(output, f + ext) for f in model_names)
    else:
        output_fnames = (output, )

    res = []
    for model_name, input_fname, output_fname in zip(model_names, input_fnames, output_fnames):
        p = Pharmacophore()
        if input_fname.endswith('.pma'):
            p.load_from_pma(input_fname)
        elif input_fname.endswith('.xyz'):
            p.load_from_xyz(input_fname)
        # skip models with less number of features with distinct coordinates that given
        if min_features is not None and len(set(xyz for label,xyz in p.get_feature_coords())) < min_features:
            continue
        p.update(bin_step=bin_step)
        fp = p.get_fp()
        res.append(Model(model_name, fp, p, output_fname))

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
    if not output_sdf:
        for mol_name, stereo_id, conf_id, out_fname in results:
            with open(out_fname, 'at') as f:
                f.write('\t'.join((mol_name, str(stereo_id), str(conf_id))) + '\n')
    else:
        for mol_name, stereo_id, conf_id, out_fname, matrix in results:
            m = db.get_mol(mol_name)[stereo_id]
            AllChem.TransformMol(m, matrix, conf_id)
            m.SetProp('_Name', f'{mol_name}-{stereo_id}-{conf_id}')
            with open(out_fname, 'a') as f:
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

    db = DB(db_fname)
    bin_step = db.get_bin_step()
    models = read_models(queries, output, output_sdf, bin_step, min_features)   # return list of Model namedtuples
    for model in models:
        if os.path.isfile(model.output_filename):
            os.remove(model.output_filename)

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
