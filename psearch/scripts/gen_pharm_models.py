#!/usr/bin/env python3
# author          : Alina Kutlushina
# date            : 10.01.2019
# license         : BSD-3
# ==============================================================================

import os.path
import sys
import time
import shelve
import argparse
import numpy as np
import pandas as pd
from collections import defaultdict
from pmapper.pharmacophore import Pharmacophore


def create_parser():
    parser = argparse.ArgumentParser(
        description='Iteratively create ligand-based pharmacophore models.')
    parser.add_argument('-p', '--project_dir', default=None,
                        help='path to a project dir. It will work if you used the standart file paths')
    parser.add_argument('-adb', '--in_active_database', metavar='active.db', required=True,
                        help='input SQL database file with active compounds')
    parser.add_argument('-idb', '--in_inactive_database', metavar='inactive.db', required=True,
                        help='input SQL database file with active compounds')
    parser.add_argument('-ats', '--in_active_trainset', metavar='active_training_set.txt', required=True,
                        help='txt file with information about active models: '
                             'model, hash, stereo, nact, ninact, nact/ninact, conf_id, feature_ids')
    parser.add_argument('-its', '--in_inactive_trainset', metavar='inactive_training_set.txt', required=True,
                        help='txt file with information about active models: '
                             'model, hash, stereo, nact, ninact, nact/ninact, conf_id, feature_ids')
    parser.add_argument('-o', '--output_path', metavar='output/path', required=False, default=None,
                        help='output path to the models of pharmacophores. '
                             'If None, the path will be generated automatically.')
    parser.add_argument('-tol', '--tolerance', default=0,
                        help='tolerance volume for the calculation of the stereo sign. If the volume of the '
                             'tetrahedron created by four points less than tolerance then those points are considered '
                             'lying on the same plane (flat; stereo sign is 0).')
    parser.add_argument('-u', '--upper', default=None,
                        help='number of features of output models')
    parser.add_argument('-l', '--lower', default=3,
                        help='number of features of input models')
    parser.add_argument('-b', '--bin_step', default=1,
                        help='binning step.')
    parser.add_argument('-mf', '--min_feature', default=None,
                        help='starting from this number of centers of pharmacophore models, models will be preserved.')
    return parser


def _keep_best_models(df, df_sub_act, df_sub_inact, save_files):
    df_sub_act = pd.merge(df_sub_act, df[['hash']], on='hash', how='inner').reset_index(drop=True)
    if not df_sub_inact.empty:
        df_sub_inact = pd.merge(df_sub_inact, df[['hash']], on='hash', how='inner').reset_index(drop=True)

    if save_files:
        path_internal = os.path.join(save_files[0], 'internal_statistics_{}_pharm{}.txt'.format(save_files[1], save_files[2]))
        path_sub_act = os.path.join(save_files[0], 'ph_active_{}_pharm{}.txt'.format(save_files[1], save_files[2]))
        df.to_csv(path_internal, index=None, sep='\t')
        df_sub_act.to_csv(path_sub_act, index=None, sep='\t')
        if not df_sub_inact.empty:
            path_sub_inact = os.path.join(save_files[0], 'ph_inactive_{}_pharm{}.txt'.format(save_files[1], save_files[2]))
            df_sub_inact.to_csv(path_sub_inact, index=None, sep='\t')
    return df_sub_act, df_sub_inact


# generator return mol_name, conf_id, hash, labels
def _gen_quadruplets(dict_ph, lower, tol):
    for mol_name in dict_ph.keys():
        for stereo_id, pharms in dict_ph[mol_name].items():
            for conf_id, pharm in enumerate(pharms):
                if pharm:
                    for hash, labels in pharm.iterate_pharm(lower, lower, tol):
                        yield mol_name, '\t'.join(map(str, [stereo_id, conf_id])), hash, labels


# generator return mol_name, conf_id, hash, labels
def _plus_one_feature(db_ph, df_sub):
    for mol_name in set(df_sub['mol_name']):
        for conf_id in set(df_sub[df_sub['mol_name'] == mol_name]['conf_id']):
            list_ids = [tuple(map(int, l.split(','))) for l in df_sub[(df_sub['mol_name'] == mol_name) &
                                                                      (df_sub['conf_id'] == conf_id)]['feature_ids']]
            sid, cid = conf_id.split()
            pharm = db_ph[mol_name][int(sid)][int(cid)]
            if pharm:
                for hash, labels in pharm.iterate_pharm1(list_ids):
                    yield mol_name, conf_id, hash, labels


# return type DataFrame: columns=['hash', 'count', 'mol_name', 'conf_id', 'feature_ids']
def gen_models(def_generator):
    dct = defaultdict(list)
    for mol_name, conf_id, hash, labels in def_generator:
        dct['mol_name'].append(mol_name)
        dct['conf_id'].append(conf_id)
        dct['hash'].append(hash)
        dct['feature_ids'].append(','.join(map(str, labels)))
    df = pd.DataFrame(dct)
    if df.empty:
        return df

    count_df = df.drop_duplicates(subset=['mol_name', 'hash'])
    count_df = count_df.groupby(['hash'], sort=True).size().reset_index(name='count')
    df = pd.merge(df, count_df, on='hash', how='right')
    df = df.sort_values(by=['count', 'hash'], ascending=False)
    return df[['hash', 'count', 'mol_name', 'conf_id', 'feature_ids']]


# return dictionary of coords: {mol_id: {stereo_id: [pharm1, pharm2, ..], stereo_id2: [..., ...]} }
def load_pharmacophores(in_db, in_training_set, bin_step):
    # read molecule names that in training set
    mol_names = [name.strip().split('\t')[1] for name in open(in_training_set).readlines()]
    conf_coords = defaultdict(dict)

    not_in_db = []
    with shelve.open(in_db, 'r') as db:
        for mol_id in mol_names:
            try:
                for sid in db[mol_id]:
                    conf_coords[mol_id] = {sid: []}
                    for coord in db[mol_id][sid]['ph']:
                        pharm = Pharmacophore(bin_step=bin_step, cached=True)
                        pharm.load_from_feature_coords(coord)
                        conf_coords[mol_id][sid].append(pharm)
            except KeyError:
                not_in_db.append(mol_id)

    if len(not_in_db) != 0:
        return '{} is/are not in {}'.format(', '.join(map(str, not_in_db)), os.path.basename(in_db)), 0
    else:
        return conf_coords, ''


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
def save_models_pma(db_ph, df_sub_act, path_pma, cluster_num, num_ids):
    df_sub_act = df_sub_act.drop_duplicates(subset=['hash'])
    i = 0
    for hash, count, mol_name, conf_id, feature_ids in df_sub_act.values:
        sid, cid = conf_id.split()
        pharm = db_ph[mol_name][int(sid)][int(cid)]
        if pharm:
            pharm.save_to_pma(os.path.join(path_pma, '{}_pharm{}_{}.pma'.format(cluster_num, num_ids, i)),
                              tuple(map(int, feature_ids.split(','))))
            i += 1
    return i


def gen_pharm_models(project_dir, in_adb, in_indb, act_trainset, inact_trainset, out_pma,
                     tolerance, lower, upper, bin_step, min_feature, save_models, save_files=False):
    time_start = time.time()
    # lower - number of model's features
    cluster_num = os.path.basename(act_trainset).split('.')[0].split('_')[1]  # number of cluster
    if cluster_num == 'centroid':
        clust_strategy = 1
    else:
        clust_strategy = 2

    db_act, not_in_db = load_pharmacophores(in_adb, act_trainset, bin_step)
    if not_in_db == 0:
        return db_act, 0

    db_inact, not_in_db = load_pharmacophores(in_indb, inact_trainset, bin_step)
    if not_in_db == 0:
        return db_inact, 0

    df_sub_act = gen_models(_gen_quadruplets(db_act, lower, tolerance))
    if db_inact:
        df_sub_inact = gen_models(_gen_quadruplets(db_inact, lower, tolerance))
    else:
        df_sub_inact = pd.DataFrame()
    df = calc_internal_stat(df_sub_act[['hash', 'count']].drop_duplicates(subset=['hash']), df_sub_inact,
                            act_trainset, clust_strategy)
    if df.empty:
        return 'no 3-points pharmacophore models', 0

    if save_files:
        path_files = os.path.join(project_dir, 'files')
        if not os.path.exists(path_files):
            os.mkdir(path_files)
        save_files = [path_files, cluster_num, lower]

    df_sub_act, df_sub_inact = _keep_best_models(df, df_sub_act, df_sub_inact, save_files)

    if save_models and lower >= min_feature:
        n = save_models_pma(db_act, df_sub_act, out_pma, cluster_num, lower)

    while True:
        if lower == upper:
            break

        lower += 1
        df_sub_act_ahead = gen_models(_plus_one_feature(db_act, df_sub_act))
        if not df_sub_inact.empty:
            df_sub_inact_ahead = gen_models(_plus_one_feature(db_inact, df_sub_inact))
        else:
            df_sub_inact_ahead = pd.DataFrame()
        df = calc_internal_stat(df_sub_act_ahead[['hash', 'count']].drop_duplicates(subset=['hash']),
                                df_sub_inact_ahead, act_trainset, clust_strategy)
        if df.empty:
            lower -= 1
            break
        if save_files:
            save_files = [os.path.join(project_dir, 'files'), cluster_num, lower]

        df_sub_act, df_sub_inact = _keep_best_models(df, df_sub_act_ahead, df_sub_inact_ahead, save_files)
        if save_models and lower >= min_feature:
            n = save_models_pma(db_act, df_sub_act, out_pma, cluster_num, lower)

    num_mdls = save_models_pma(db_act, df_sub_act, out_pma, cluster_num, lower)

    sys.stderr.write('{}: extract {} models passed ({}s)\n'.format(cluster_num, num_mdls, round(time.time() - time_start, 3)))
    return out_pma, lower


def entry_point():
    parser = create_parser()
    args = vars(parser.parse_args())
    for o, v in args.items():
        if o == "project_dir": project_dir = v
        if o == "in_active_database": adb = v
        if o == "in_inactive_database": indb = v
        if o == "in_active_trainset": act_trainset = v
        if o == "in_inactive_trainset": inact_trainset = v
        if o == "output_path": out_pma = v
        if o == "tolerance": tolerance = float(v)
        if o == "lower": lower = int(v)
        if o == "upper": upper = int(v) if v is not None else None
        if o == "bin_step": bin_step = int(v)
        if o == "min_feature": min_feature = int(v)

    if not out_pma:
        out_pma = os.path.join(project_dir, 'models')

    if min_feature:
        save_models = True
    else:
        save_models = False

    gen_pharm_models(project_dir=project_dir,
                     in_adb=adb,
                     in_indb=indb,
                     act_trainset=act_trainset,
                     inact_trainset=inact_trainset,
                     out_pma=out_pma,
                     tolerance=tolerance,
                     lower=lower,
                     upper=upper,
                     bin_step=bin_step,
                     min_feature=min_feature,
                     save_models=save_models)


if __name__ == '__main__':
    entry_point()
