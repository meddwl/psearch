#!/usr/bin/env python3
# author          : Pavel Polishchuk
# date            : 23.08.2019
# license         : BSD-3
# ==============================================================================

import os
import argparse
import marshal
import sqlite3
from collections import defaultdict, namedtuple
from pmapper.pharmacophore import Pharmacophore
from multiprocessing import Pool, Lock
from functools import partial

Model = namedtuple('Model', ['name', 'fp', 'pharmacophore', 'output_filename'])
Conformer = namedtuple('Conformer', ['id', 'fp', 'pharmacophore'])


def load_confs(mol_name, db_fname):
    connection = sqlite3.connect(db_fname)
    cur = connection.cursor()

    # load all fp for conformers of a molecule
    cur.execute("SELECT conf_id, fp FROM conformers WHERE mol_name = ?", (mol_name,))
    data = cur.fetchall()   # (conf_id, fp)
    data = {conf_id: marshal.loads(fp) for conf_id, fp in data}   # {conf_id: fp_unpacked}

    # load feature coordinates of all conformers
    sql = "SELECT conf_id, feature_label, x, y, z FROM feature_coords WHERE conf_id IN (%s)" % ','.join(map(str, data.keys()))
    cur.execute(sql)
    res = cur.fetchall()    # (conf_id, feature_label, x, y, z)
    feature_coords = defaultdict(list)
    for r in res:
        feature_coords[r[0]].append((r[1], tuple(r[2:])))   # {conf_id: [(label, (x, y, z)), (label, (x, y, z)), ...]}

    # get bin step from DB
    cur.execute("SELECT bin_step FROM settings")
    bin_step = cur.fetchone()[0]

    res = []
    for conf_id in data.keys():
        p = Pharmacophore(bin_step=bin_step)
        p.load_from_feature_coords(feature_coords[conf_id])
        res.append(Conformer(conf_id, data[conf_id], p))

    return res   # list of Conformers


def get_comp_names_from_db(db_fname):
    conn = sqlite3.connect(db_fname)
    cur = conn.cursor()
    cur.execute("SELECT DISTINCT(mol_name) FROM conformers")
    mol_names = tuple(i[0] for i in cur.fetchall())
    return mol_names


def read_models(queries, output, is_output_sdf):

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
        fp = p.get_fp()
        res.append(Model(model_name, fp, p, output_fname))

    return res


def screen(mol_name, db_conn, models, input_sdf, match_first_conf):

    def compare_fp(query_fp, fp):
        return (query_fp & fp) == query_fp

    get_transform_matrix = input_sdf is not None

    confs = load_confs(mol_name, db_conn)

    output = []
    for model in models:
        for conf in confs:
            if compare_fp(model.fp, conf.fp):
                res = conf.pharmacophore.fit_model(model.pharmacophore, get_transform_matrix=get_transform_matrix)
                if res:
                    if get_transform_matrix:
                        output.append((mol_name, conf.id, model.output_filename, res[1]))
                    else:
                        output.append((mol_name, conf.id, model.output_filename))
                    if match_first_conf:
                        break
    return output


def screen_db(db_fname, queries, output, input_sdf, match_first_conf, ncpu):

    if output.endswith('.txt') or output.endswith('.sdf'):
        if not os.path.exists(os.path.dirname(output)):
            os.makedirs(os.path.dirname(output), exist_ok=True)
    else:
        if not os.path.exists(output):
            os.makedirs(output, exist_ok=True)

    is_sdf_output = input_sdf is not None
    models = read_models(queries, output, is_sdf_output)   # return list of Model namedtuples
    for model in models:
        if os.path.isfile(model.output_filename):
            os.remove(model.output_filename)

    comp_names = get_comp_names_from_db(db_fname)
    if ncpu == 1:
        for res in screen(mol_name=comp_names, db_conn=db_fname, models=models, input_sdf=input_sdf, match_first_conf=match_first_conf):
            if not is_sdf_output:
                for mol_name, conf_id, out_fname in res:
                    with open(out_fname, 'at') as f:
                        f.write('{}\t{}\n'.format(mol_name, conf_id))
    else:
        p = Pool(ncpu)
        for res in p.imap_unordered(partial(screen, db_conn=db_fname, models=models, input_sdf=input_sdf, match_first_conf=match_first_conf), comp_names, chunksize=10):
            if not is_sdf_output:
                for mol_name, conf_id, out_fname in res:
                    with open(out_fname, 'at') as f:
                        f.write('{}\t{}\n'.format(mol_name, conf_id))
        p.close()


def entry_point():
    parser = argparse.ArgumentParser(description='Screen SQLite DB with compounds against pharmacophore queries.',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-d', '--dbname', metavar='input_active.db, input_inactive.db',
                        type=str, required=True,
                        help='input SQLite database with generated conformers.')
    parser.add_argument('-q', '--query', metavar='model.pma', required=True, type=str, nargs='+',
                        help='pharmacophore model or models or a directory path. If a directory is specified all '
                             'pma- and xyz-files will be used for screening as pharmacophore models.')
    parser.add_argument('-o', '--output', required=True, type=str,
                        help='path to an output text (.txt) file which will store names of compounds fit the model(s). '
                             'If input_sdf argument is specified the output should be an sdf file (.sdf) to store '
                             'conformers fitted to a model. In the case multiple models were supplied this '
                             'should be path to a directory where output files will be created to store '
                             'screening results. Type of the output is recognized by file extension. '
                             'Existed output files will be overwritten.')
    parser.add_argument('--input_sdf', metavar='input.sdf', default=None, type=str,
                        help='sdf file with conformers used for creation of SQLite DB. Should be specified if '
                             'conformers fitted to model should be returned.')
    parser.add_argument('--conf', action='store_true', default=False, type=bool,
                        help='return all conformers matches as separate hits in a hit list.')
    parser.add_argument('-c', '--ncpu', metavar='INTEGER', default=1, type=int,
                        help='number of cores to use.')

    args = parser.parse_args()

    screen_db(db_fname=args.dbname,
              queries=args.query,
              output=args.output,
              input_sdf=args.input_sdf,
              match_first_conf=not args.conf,
              ncpu=args.ncpu)


if __name__ == '__main__':
    entry_point()
