#!/usr/bin/env python3
# author          : Alina Kutlushina
# date            : 01.05.2018
# license         : BSD-3
#==============================================================================

import os
import sys
import time
import argparse
import pandas as pd
import sqlite3 as sqlite
from collections import defaultdict
from pmapper.pharmacophore import Pharmacophore


def load_pharmacophores(in_db, in_training_set):
#     start_time_r = time.time()
    mol_names = [name.strip().split('\t')[0] for name in open(in_training_set).readlines()]
 
    confs_pharm = defaultdict(dict)
    with sqlite.connect(in_db) as con:
        for mol_name in mol_names:
            cur = con.cursor()
            cur.execute("SELECT conf_id, feature_label, x, y, z FROM feature_coords WHERE conf_id IN "
                        "(SELECT conf_id from conformers WHERE mol_name = ?)", (mol_name,))
            res = cur.fetchall()
            confs = defaultdict(list)
            for r in res:
                confs[r[0]].append((r[1], tuple(r[2:])))  # dict(conf_id: (feature_label, x, y, z))
            for conf_id, coord in confs.items():
                confs_pharm[mol_name][conf_id] = Pharmacophore()
                confs_pharm[mol_name][conf_id].load_from_feature_coords(coord)
#     sys.stderr.write('(%is)\n' % (time.time() - start_time_r))
    return confs_pharm


def get_data(confs_pharm, p_num, tol):
    for conf_id, pharm in confs_pharm.items():
        yield conf_id, pharm, p_num, tol


def gen_subpharm(values):
    conf_id, pharm, p_num, tol = values
    subph_inf = []
    for hash, labels in pharm.iterate_pharm(p_num, p_num, tol):
        subph_inf.append((conf_id, hash, labels))
    return subph_inf


def main(in_db, fname, tol, p_num, freq, out_fname=None):

    pharms = load_pharmacophores(os.path.abspath(in_db), fname)

    if out_fname is None:
        out_fname = os.path.join(os.path.dirname(fname),
                             'ph_%s_pharm%i.txt' % (os.path.splitext(os.path.basename(fname))[0], p_num))

    start_time = time.time()
    dct = defaultdict(list)
    for mol_name, confs_pharm in pharms.items():
        for conf_id, pharm in confs_pharm.items():
            for hash, labels in pharm.iterate_pharm(p_num, p_num, tol):
                dct['hash'].append(hash)
                dct['mol_name'].append(mol_name)
                dct['conf_id'].append(conf_id)
                dct['feature_ids'].append(','.join(map(str, labels)))

    df = pd.DataFrame(dct)
    df = df.drop_duplicates(subset=['mol_name', 'hash'])
    count_df = df.groupby(['hash'], sort=True).size().reset_index(name='count')
    full_df = pd.merge(df, count_df, on='hash', how='right')
    full_df = full_df.sort_values(by=['count', 'hash'], ascending=False)
    full_df = full_df[['hash', 'count', 'mol_name', 'conf_id', 'feature_ids']]
    full_df.to_csv(out_fname, index=None, sep='\t')

    sys.stderr.write('%s: features %i passed (%is)\n' % (os.path.basename(fname), p_num, time.time() - start_time))
    sys.stderr.flush()

    return out_fname


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-d', '--input_db', metavar='input.db', required=True,
                        help='SQLite DB with pharmacophores (feature coordinates).')
    parser.add_argument('-ts', '--file_trainset', default=None,
                        help='txt file with name of molecules for training set')
    parser.add_argument('-t', '--tolerance', metavar='VALUE', default=0,
                        help='tolerance used for calculation of a stereoconfiguration sign.')
    parser.add_argument('-num', '--num_of_reatures', metavar='4', type=int, default=0,
                        help='number of features used for generation of subpharmacophores.')
    parser.add_argument('-f', '--freq', metavar='0', default=None,
                        help='minimal frequency of pharmacophore occurrence in input compounds to keep them '
                             'in output file. Default: None (pharmacophores occurred in at least a half of '
                             'input compounds will be stored).')

    args = vars(parser.parse_args())
    for o, v in args.items():
        if o == "input_db": in_db = v
        if o == "file_trainset": file_ts = v
        if o == "tolerance": tol = v
        if o == "num_of_reatures": num = int(v)
        if o == "freq": freq = int(v) if v else None

    main(in_db=in_db,
         fname=file_ts,
         tol=tol,
         p_num=num,
         freq=freq)
