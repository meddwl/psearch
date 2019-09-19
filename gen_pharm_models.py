#!/usr/bin/env python3
# author          : Alina Kutlushina
# date            : 10.01.2019
# license         : BSD-3
# ==============================================================================

import os.path
import sys
import time
import argparse
import numpy as np
import pandas as pd
import sqlite3 as sql
from collections import defaultdict
from pharmacophore import Pharmacophore


def _make_dir(new_path):
    for path in [os.path.split(new_path)[0], new_path]:
        if not os.path.exists(path):
            os.mkdir(path)


def _keep_best_models(df, df_sub_act, df_sub_inact, df_ph_act, df_ph_inact, save_files):
    df_sub_act = pd.merge(df_sub_act, df[['hash']], on='hash', how='inner').reset_index(drop=True)
    df_ph_act = pd.merge(df_ph_act, df_sub_act[['conf_id']].drop_duplicates(subset=['conf_id']), on='conf_id',
                         how='inner').reset_index(drop=True)
    if not df_ph_inact.empty:
        df_sub_inact = pd.merge(df_sub_inact, df[['hash']], on='hash', how='inner').reset_index(drop=True)
        df_ph_inact = pd.merge(df_ph_inact, df_sub_inact[['conf_id']].drop_duplicates(subset=['conf_id']),
                           on='conf_id', how='inner').reset_index(drop=True)
    if save_files:
        path_internal = os.path.join(save_files[0], 'internal_statistics_{}_pharm{}.txt'.format(save_files[1], save_files[2]))
        path_sub_act = os.path.join(save_files[0], 'ph_active_{}_pharm{}.txt'.format(save_files[1], save_files[2]))
        df.to_csv(path_internal, index=None, sep='\t')
        df_sub_act.to_csv(path_sub_act, index=None, sep='\t')
        if not df_sub_inact.empty:
            path_sub_inact = os.path.join(save_files[0], 'ph_inactive_{}_pharm{}.txt'.format(save_files[1], save_files[2]))
            df_sub_inact.to_csv(path_sub_inact, index=None, sep='\t')
    return df_sub_act, df_sub_inact, df_ph_act, df_ph_inact


# generator return mol_name, conf_id, hash, labels
def _gen_quadruplets(df_ph, lower, tol):
    for mol_name, conf_id, pharm in zip(df_ph['mol_name'], df_ph['conf_id'], df_ph['pharm']):
        if pharm:
            for hash, labels in pharm.iterate_pharm(lower, lower, tol):
                yield mol_name, conf_id, hash, labels


# generator return mol_name, conf_id, hash, labels
def _plus_one_feature(df_ph, df_sub):
    for mol_name, conf_id, pharm in zip(df_ph['mol_name'], df_ph['conf_id'], df_ph['pharm']):
        list_ids = df_sub[df_sub['conf_id'] == conf_id]
        list_ids = [tuple(map(int, l.split(','))) for l in list_ids['feature_ids']]
        if pharm:
            for hash, labels in pharm.iterate_pharm1(list_ids):
                yield mol_name, conf_id, hash, labels


# return type DataFrame: columns=['hash', 'count', 'mol_name', 'conf_id', 'feature_ids']
def gen_models(def_generator, df_0):
    dct = defaultdict(list)
    for mol_name, conf_id, hash, labels in def_generator:
        dct['hash'].append(hash)
        dct['mol_name'].append(mol_name)
        dct['conf_id'].append(conf_id)
        dct['feature_ids'].append(','.join(map(str, labels)))
    df = pd.DataFrame(dct)
    if df.empty:
        print(df_0[:2])
        return df_0, df
    count_df = df.drop_duplicates(subset=['mol_name', 'hash'])
    count_df = count_df.groupby(['hash'], sort=True).size().reset_index(name='count')
    df = pd.merge(df, count_df, on='hash', how='right')
    df = df.sort_values(by=['count', 'hash'], ascending=False)
    return df_0, df[['hash', 'count', 'mol_name', 'conf_id', 'feature_ids']]


# return DataFrame of pharmacophore representation molecules: columns=['mol_name', 'conf_id', 'pharm']
def load_pharmacophores(in_db, in_training_set):
    mol_names = [name.strip().split('\t')[1] for name in open(in_training_set).readlines()]
    confs_pharm = defaultdict(list)
    with sql.connect(in_db) as con:
        cur = con.cursor()
        cur.execute("SELECT bin_step FROM settings")
        db_bin_step = cur.fetchone()[0]
        for mol_name in mol_names:
            cur.execute("SELECT conf_id, feature_label, x, y, z FROM feature_coords WHERE conf_id IN "
                        "(SELECT conf_id from conformers WHERE mol_name = ?)", (mol_name,))
            res = cur.fetchall()
            confs = defaultdict(list)
            for r in res:
                confs[r[0]].append((r[1], tuple(r[2:])))  # dict(conf_id: (feature_label, x, y, z))
            for conf_id, coord in confs.items():
                p = Pharmacophore(bin_step=db_bin_step, cached=True)
                p.load_from_feature_coords(coord)
                confs_pharm['mol_name'].append(mol_name)
                confs_pharm['conf_id'].append(conf_id)
                confs_pharm['pharm'].append(p)
    return pd.DataFrame(confs_pharm)


# return type DataFrame
def strategy_extract_trainset(df, clust_strategy):
    if clust_strategy == 2:
        df = df.sort_values(by=['recall', 'F2', 'F05'], ascending=False)
        df = df.reset_index(drop=True)
        if df['F2'].iloc[0] == 1.0:
            df = df[(df['recall'] == 1.0) & (df['F2'] == 1.0)]
        elif df[df['F2'] >= 0.8].shape[0] <= 100:
            df = df[(df['recall'] == 1) & (df['F2'] >= 0.8)]
        else:
            df = df[(df['recall'] == 1) & (df['F2'] >= df['F2'].loc[100])]

    elif clust_strategy == 1:
        df = df.sort_values(by=['recall', 'F05', 'F2'], ascending=False)
        df = df.reset_index(drop=True)
        if df[df['F05'] >= 0.8].shape[0] <= 100:
            df = df[df['F05'] >= 0.8]
        else:
            df = df[df['F05'] >= df['F05'].loc[100]]
    return df


# return type DataFrame: columns=['hash', 'TP', 'FP', 'precision', 'recall', 'F2', 'F05']
def calc_internal_stat(df_act, df_inact, act_trainset, clust_strategy):
    df_act = df_act.rename(columns={'count': 'TP'})

    if not df_inact.empty:
        df_inact = df_inact[['hash', 'count']].drop_duplicates(subset=['hash'])
        df_inact = df_inact.rename(columns={'count': 'FP'})
        df = pd.merge(df_act, df_inact, on=['hash'], how='left')
        df.loc[df['FP'].isnull(), 'FP'] = 0
    else:
        df = df_act
        df['FP'] = [0 for x in range(df.shape[0])]

    df['precision'] = round(df['TP'] / (df['TP'] + df['FP']), 3)
    df['recall'] = round(df['TP'] / (len(open(act_trainset).readlines())), 3)
    df['F2'] = round(5 * ((df['precision'] * df['recall']) / (4 * df['precision'] + df['recall'])), 3)
    df['F05'] = round(1.25 * ((df['precision'] * df['recall']) / (0.25 * df['precision'] + df['recall'])), 3)
    df['FP'] = df['FP'].astype(np.int64)
    df = df.sort_values(by=['recall', 'F2', 'F05'], ascending=False)
    df = df[['hash', 'TP', 'FP', 'precision', 'recall', 'F2', 'F05']]
    # difference ways to check out the best models
    df = strategy_extract_trainset(df, clust_strategy)
    return df


# return None
def save_models_pma(df_ph, df_sub_act, path_pma, cluster_num, num_ids):
    time_start = time.time()
    df_sub_act = df_sub_act.drop_duplicates(subset=['hash'])
    df_ph = pd.merge(df_sub_act, df_ph, on=['mol_name', 'conf_id'], how='inner')
    i = 0
    for line in range(df_ph.shape[0]):
        pharm, ids = df_ph['pharm'].iloc[line], df_ph['feature_ids'].iloc[line]
        if pharm:
            pharm.save_to_pma(os.path.join(path_pma, '{}_pharm{}_{}.pma'.format(cluster_num, num_ids, i)),
                              tuple(map(int, ids.split(','))))
            i += 1
    sys.stderr.write('{}: extract {} models passed ({}s)\n\n'.format(
        os.path.split(path_pma)[1], i, round(time.time() - time_start, 3)))


def main(in_adb, in_indb, act_trainset, inact_trainset, out_pma, tolerance, lower, save_files=False):
    # lower - number of model's features
    cluster_num = os.path.basename(act_trainset).split('.')[0].split('_')[1]  # number of cluster
    if cluster_num == 'centroid':
        clust_strategy = 1
    else:
        clust_strategy = 2

    df_ph_act = load_pharmacophores(in_adb, act_trainset)
    df_ph_inact = load_pharmacophores(in_indb, inact_trainset)

    if df_ph_act.empty:
        return 'Molecules from the training set are not in the database.', 0

    df_ph_act, df_sub_act = gen_models(_gen_quadruplets(df_ph_act, lower, tolerance), df_ph_act)
    if not df_ph_inact.empty:
        df_ph_inact, df_sub_inact = gen_models(_gen_quadruplets(df_ph_inact, lower, tolerance), df_ph_inact)
    else:
        df_sub_inact = pd.DataFrame()
    df = calc_internal_stat(df_sub_act[['hash', 'count']].drop_duplicates(subset=['hash']),
                           df_sub_inact,
                           act_trainset, clust_strategy)
    if df.empty:
        return 'no {}-points pharmacophore models above thresholds'.format(lower), 0

    if save_files:
        path_files = os.path.join(os.path.split(os.path.dirname(in_adb))[0], 'files')
        if not os.path.exists(path_files):
            os.mkdir(path_files)
        save_files = [path_files, cluster_num, lower]

    df_sub_act, df_sub_inact, df_ph_act, df_ph_inact = _keep_best_models(
        df, df_sub_act, df_sub_inact, df_ph_act, df_ph_inact, save_files)

    while True:
        lower += 1

        print(lower, df_sub_act[:3], sep='\n')

        df_sub_act_0, df_sub_act = gen_models(_plus_one_feature(df_ph_act, df_sub_act[['conf_id', 'feature_ids']]), df_sub_act)
        if not df_sub_inact.empty:
            df_sub_inact_0, df_sub_inact = gen_models(_plus_one_feature(df_ph_inact, df_sub_inact[['conf_id', 'feature_ids']]), df_sub_inact)
        df = calc_internal_stat(df_sub_act[['hash', 'count']].drop_duplicates(subset=['hash']),
                                df_sub_inact,
                                act_trainset, clust_strategy)
        if df.empty:
            lower -= 1
            break

        if save_files:
            save_files = [os.path.join(os.path.split(os.path.dirname(in_adb))[0], 'files'), cluster_num, lower]

        df_sub_act, df_sub_inact, df_ph_act, df_ph_inact = _keep_best_models(
            df, df_sub_act, df_sub_inact, df_ph_act, df_ph_inact, save_files)

    path_pma = os.path.join(out_pma, '{}_ph{}'.format(cluster_num, lower))
    _make_dir(path_pma)
    save_models_pma(df_ph_act, df_sub_act_0, path_pma, cluster_num, lower)
    return path_pma, lower


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='ilter out duplicated hash-stereo pair for each compound ID (without stereo)')
    parser.add_argument('-adb', '--in_active_database', metavar='active.db', required=True,
                        help='input SQL database file with active compounds')
    parser.add_argument('-idb', '--in_inactive_database', metavar='inactive.db', required=True,
                        help='input SQL database file with active compounds')
    parser.add_argument('-ats', '--in_active_trainset', metavar='active_training_set.txt', required=True,
                        help='txt file with information adout active models: '
                             'model, hash, stereo, nact, ninact, nact/ninact, conf_id, feature_ids')
    parser.add_argument('-its', '--in_inactive_trainset', metavar='inactive_training_set.txt', required=True,
                        help='txt file with information adout active models: '
                             'model, hash, stereo, nact, ninact, nact/ninact, conf_id, feature_ids')
    parser.add_argument('-o', '--output_path', metavar='output/path', required=False, default=None,
                        help='output path to the models of pharmacophores. '
                             'If None, the path will be generated automatically.')
    parser.add_argument('-tol', '--tolerance', default=0,
                        help='tolerance volume for the calculation of the stereo sign. If the volume of the '
                             'tetrahedron created by four points less than tolerance then those points are considered '
                             'lying on the same plane (flat; stereo sign is 0).')
    parser.add_argument('-l', '--lower', default=4,
                        help='number of features of input models')

    args = vars(parser.parse_args())
    for o, v in args.items():
        if o == "in_active_database": adb = v
        if o == "in_inactive_database": indb = v
        if o == "in_active_trainset": act_trainset = v
        if o == "in_inactive_trainset": inact_trainset = v
        if o == "output_path": out_pma = v
        if o == "tolerance": tolerance = float(v)
        if o == "lower": lower = int(v)

    if not out_pma:
        out_pma = os.path.join(os.path.split(os.path.dirname(adb))[0], 'models')

    main(in_adb=adb,
         in_indb=indb,
         act_trainset=act_trainset,
         inact_trainset=inact_trainset,
         out_pma=out_pma,
         tolerance=tolerance,
         lower=lower)
