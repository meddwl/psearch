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
from psearch.database import DB
from pmapper.pharmacophore import Pharmacophore


def create_parser():
    parser = argparse.ArgumentParser(description='Iteratively create ligand-based pharmacophore models.',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-db', '--path_database', metavar='active.db', required=True,
                        help='input SQL database file with active compounds')
    parser.add_argument('-ats', '--active_trainset', metavar='active_training_set.txt', required=True,
                        help='txt file with information about active models: '
                             'model, hash, stereo, nact, ninact, nact/ninact, conf_id, feature_ids')
    parser.add_argument('-its', '--inactive_trainset', metavar='inactive_training_set.txt', required=True,
                        help='txt file with information about active models: '
                             'model, hash, stereo, nact, ninact, nact/ninact, conf_id, feature_ids')
    parser.add_argument('-o', '--path_models', metavar='output/path', required=False, default=None,
                        help='output path to the models of pharmacophores. '
                             'If None, the path will be generated automatically.')
    parser.add_argument('-tol', '--tolerance', type=int,  default=0,
                        help='tolerance volume for the calculation of the stereo sign. If the volume of the '
                             'tetrahedron created by four points less than tolerance then those points are considered '
                             'lying on the same plane (flat; stereo sign is 0).')
    parser.add_argument('-u', '--upper', type=int,  default=None,
                        help='number of features of output models')
    parser.add_argument('-l', '--lower', type=int, default=3,
                        help='number of features of input models')
    parser.add_argument('-b', '--bin_step', type=int, default=1,
                        help='binning step.')
    parser.add_argument('-mf', '--min_feature', type=int, default=0,
                        help='starting from this number of centers of pharmacophore models, models will be preserved.')
    return parser


def _keep_best_models(df, df_pharm_act, df_pharm_inact, df_sub_act, df_sub_inact, save_files):

    df_sub_act = df_sub_act[df_sub_act['hash'].isin(set(df['hash']))].reset_index(drop=True)
    df_pharm_act = df_pharm_act[df_pharm_act['conf_id'].isin(set(df_sub_act['conf_id']))].reset_index(drop=True)
    if not df_sub_inact.empty:
        df_sub_inact = df_sub_inact[df_sub_inact['hash'].isin(set(df['hash']))].reset_index(drop=True)
        df_pharm_inact = df_pharm_inact[df_pharm_inact['conf_id'].isin(set(df_sub_act['conf_id']))].reset_index(drop=True)

    if save_files:
        path_internal = os.path.join(save_files[0], 'internal_statistics_{}_f{}.txt'.format(save_files[1], save_files[2]))
        path_sub_act = os.path.join(save_files[0], 'ph_active_{}_f{}.txt'.format(save_files[1], save_files[2]))
        df.to_csv(path_internal, index=None, sep='\t')
        df_sub_act.to_csv(path_sub_act, index=None, sep='\t')
        if not df_sub_inact.empty:
            path_sub_inact = os.path.join(save_files[0], 'ph_inactive_{}_f{}.txt'.format(save_files[1], save_files[2]))
            df_sub_inact.to_csv(path_sub_inact, index=None, sep='\t')
    return df_pharm_act, df_pharm_inact, df_sub_act, df_sub_inact


# generator return mol_name, conf_id, hash, labels
def _gen_quadruplets(db, pp_train_set, lower, tol, bin_step):
    mol_names = [name.strip().split()[1] for name in open(pp_train_set).readlines()]
    for mol_name in mol_names:
        try:
            dict_coords = db.get_pharm(mol_name)
        except KeyError:
            sys.exit(f"Unexpect molecule name. This {mol_name} is not in the input databased")
        for isomer_id, list_coords in dict_coords.items():
            for conf_id, coord in enumerate(list_coords):
                pharm = Pharmacophore(bin_step=bin_step, cached=True)
                pharm.load_from_feature_coords(coord)
                if pharm:
                    conf_id = f'{isomer_id}_{conf_id}'
                    for hash, labels in pharm.iterate_pharm(lower, lower, tol):
                        yield mol_name, conf_id, hash, labels


# generator return mol_name, conf_id, hash, labels
def _plus_one_feature(db, df_pharm, df_sub, bin_step):
    for mol_name, conf_id in df_pharm.values:
        feature_ids = set(df_sub[(df_sub['mol_name'] == mol_name) & (df_sub['conf_id'] == conf_id)]['feature_ids'])
        label_ids = [tuple(map(int, lbls.split(','))) for lbls in feature_ids]
        sid, cid = conf_id.split('_')
        pharm = Pharmacophore(bin_step=bin_step, cached=True)
        pharm.load_from_feature_coords(db.get_pharm(mol_name)[int(sid)][int(cid)])
        if pharm:
            for hash, labels in pharm.iterate_pharm1(label_ids):
                yield mol_name, conf_id, hash, labels


# return type DataFrame: columns=['hash', 'count', 'mol_name', 'conf_id', 'feature_ids']
def gen_models(def_generator):
    data = []
    for list_data in def_generator:
        data.append(list(list_data[:3]) + [','.join(map(str, list_data[3]))])
    df = pd.DataFrame(data, columns=['mol_name', 'conf_id', 'hash', 'feature_ids'])
    if df.empty:
        return df
    count_df = df.drop_duplicates(subset=['mol_name', 'hash'])
    count_df = count_df.groupby(['hash'], sort=True).size().reset_index(name='count')
    df = pd.merge(df, count_df, on='hash', how='right')
    df = df.sort_values(by=['count', 'hash'], ascending=False)
    return df[['hash', 'count', 'mol_name', 'conf_id', 'feature_ids']]


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
    if not df_inact.empty:
        df = pd.merge(df_act, df_inact, on=['hash'], how='left')
        df.loc[df['FP'].isnull(), 'FP'] = 0
    else:
        df = df_act
        df['FP'] = [0] * df.shape[0]

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
def save_models_pma(db, df_sub_act, path_pma, bin_step, cluster_num, num_ids):
    data = df_sub_act.drop_duplicates(subset=['hash']).values
    for num, (hash, count, mol_name, conf_id, feature_ids) in enumerate(data):
        isomer_id, conf_id = conf_id.split('_')
        pharm = Pharmacophore(bin_step=bin_step, cached=True)
        pharm.load_from_feature_coords(db.get_pharm(mol_name)[int(isomer_id)][int(conf_id)])
        pharm.save_to_pma(os.path.join(path_pma, f'{cluster_num}_f{num_ids}_p{num}.pma'),
                          tuple(map(int, feature_ids.split(','))))
    return len(data)


def gen_pharm_models(project_dir, in_db, act_trainset, inact_trainset, out_pma,
                     tolerance, lower, upper, bin_step, save_models, save_files=False):

    time_start = time.time()
    # lower - number of model's features
    cluster_num = os.path.splitext(os.path.basename(act_trainset))[0].split('_')[1]  # number of cluster
    if cluster_num == 'centroid':
        clust_strategy = 1
    else:
        clust_strategy = 2

    db = DB(in_db)
    df_sub_act = gen_models(_gen_quadruplets(db, act_trainset, lower, tolerance, bin_step))
    df_sub_inact = gen_models(_gen_quadruplets(db, inact_trainset, lower, tolerance, bin_step))
    df_pharm_act = df_sub_act[['mol_name', 'conf_id']].drop_duplicates(subset=['mol_name', 'conf_id'])
    df_pharm_inact = df_sub_inact[['mol_name', 'conf_id']].drop_duplicates(subset=['mol_name', 'conf_id'])
    df = calc_internal_stat(df_sub_act[['hash', 'count']].drop_duplicates(subset=['hash']).rename(columns={'count': 'TP'}),
                            df_sub_inact[['hash', 'count']].drop_duplicates(subset=['hash']).rename(columns={'count': 'FP'}),
                            act_trainset, clust_strategy)
    if df.empty:
        sys.stderr.write(f'no 3-points pharmacophore models for {cluster_num} training set\n')

    if save_files:
        path_files = os.path.join(project_dir, 'files')
        if not os.path.exists(path_files):
            os.mkdir(path_files)
        save_files = [path_files, cluster_num, lower]

    df_pharm_act, df_pharm_inact, df_sub_act, df_sub_inact = _keep_best_models(df, df_pharm_act, df_pharm_inact,
                                                                               df_sub_act, df_sub_inact, save_files)

    if save_models != 0 and lower >= save_models:
        _ = save_models_pma(db, df_sub_act, out_pma, bin_step, cluster_num, lower)

    while True:
        if lower == upper:
            break

        lower += 1
        df_sub_act_ahead = gen_models(_plus_one_feature(db, df_pharm_act, df_sub_act, bin_step))
        if df_sub_inact.empty:
            df_sub_inact_ahead = pd.DataFrame(columns=['hash', 'count', 'mol_name', 'conf_id', 'feature_ids'])
        else:
            df_sub_inact_ahead = gen_models(_plus_one_feature(db, df_pharm_inact, df_sub_inact, bin_step))

        df = calc_internal_stat(df_sub_act_ahead[['hash', 'count']].drop_duplicates(subset='hash').rename(columns={'count': 'TP'}),
                                df_sub_inact_ahead[['hash', 'count']].drop_duplicates(subset='hash').rename(columns={'count': 'FP'}),
                                act_trainset, clust_strategy)
        if df.empty:
            break
        if save_files:
            save_files = [os.path.join(project_dir, 'files'), cluster_num, lower]

        df_pharm_act, df_pharm_inact, df_sub_act, df_sub_inact = _keep_best_models(df, df_pharm_act, df_pharm_inact,
                                                                                   df_sub_act_ahead, df_sub_inact_ahead,
                                                                                   save_files)
        if save_models != 0 and lower >= save_models:
            _ = save_models_pma(db, df_sub_act, out_pma, bin_step, cluster_num, lower)

    num_models = save_models_pma(db, df_sub_act, out_pma, bin_step, cluster_num, lower - 1)
    sys.stderr.write(f'{cluster_num}: extract {num_models} models passed ({round(time.time()-time_start, 3)}s)\n')
    sys.stderr.flush()


if __name__ == '__main__':
    parser = create_parser()
    args = parser.parse_args()

    if not args.path_models:
        out_model = os.path.join(args.project_dir, 'models')
    else:
        out_model = args.path_models
    if not os.path.exists(out_model):
        os.makedirs(out_model)

    gen_pharm_models(project_dir=args.project_dir,
                     in_db=args.active_database,
                     act_trainset=args.active_trainset,
                     inact_trainset=args.inactive_trainset,
                     out_pma=out_model,
                     tolerance=args.tolerance,
                     lower=args.lower,
                     upper=args.upper if args.upper is not None else None,
                     bin_step=args.bin_step,
                     save_models=args.save_models)
