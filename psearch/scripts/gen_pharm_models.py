#!/usr/bin/env python3
# author          : Alina Kutlushina
# date            : 10.01.2019
# license         : BSD-3
# ==============================================================================

import os.path
import sys
import time
import argparse
import pandas as pd
from psearch.database import DB
from pmapper.pharmacophore import Pharmacophore


def _keep_best_models(df, df_sub, save_files, nfeatures):
    df_sub = df_sub[df_sub['hash'].isin(set(df['hash']))].reset_index(drop=True)
    if save_files:
        df.to_csv(os.path.join(save_files[0], f'internal_statistics-{save_files[1]}-f{nfeatures}.txt'), index=None, sep='\t')
        df_sub.to_csv(os.path.join(save_files[0], f'pharmacophore_models-{save_files[1]}-f{nfeatures}.txt'), index=None, sep='\t')
    return df_sub


def _gen_quadruplets(db, pp_train_set, lower, tol, bin_step):
    train_set_list = [name.strip().split() for name in open(pp_train_set).readlines()]
    for _, mol_name, activity in train_set_list:
        try:
            dict_coords = db.get_pharm(mol_name)
        except KeyError:
            sys.exit(f"Unexpect molecule name. This {mol_name} is not in the input databased")
        for isomer_id, list_coords in dict_coords.items():
            for conf_id, coord in enumerate(list_coords):
                pharm = Pharmacophore(bin_step=bin_step, cached=True)
                pharm.load_from_feature_coords(coord)
                if pharm:
                    for hash, labels in pharm.iterate_pharm(lower, lower, tol):
                        yield activity, mol_name, isomer_id, conf_id, hash, labels


def _plus_one_feature(db, df_sub, bin_step):
    df_sub_group = df_sub.groupby(['activity', 'mol_name', 'isomer_id', 'conf_id'])
    cols = ['activity', 'mol_name', 'isomer_id', 'conf_id']
    for df_group in df_sub_group:
        activity, mol_name, isomer_id, conf_id = df_group[1][cols].iloc[0]
        label_ids = [tuple(map(int, lbls.split(','))) for lbls in df_group[1]['feature_ids'].tolist()]
        pharm = Pharmacophore(bin_step=bin_step, cached=True)
        pharm.load_from_feature_coords(db.get_pharm(mol_name)[isomer_id][conf_id])
        if pharm:
            for hash, labels in pharm.iterate_pharm1(label_ids):
                yield activity, mol_name, isomer_id, conf_id, hash, labels


def gen_models(def_generator):
    data = []
    for activity, mol_name, isomer_id, conf_id, hash, labels in def_generator:
        data.append([activity, mol_name, isomer_id, conf_id, hash, ','.join(map(str, labels))])
    df = pd.DataFrame(data, columns=['activity', 'mol_name', 'isomer_id', 'conf_id', 'hash', 'feature_ids'])
    if df.empty:
        return df
    count_df = df.drop_duplicates(subset=['activity', 'mol_name', 'hash'])
    count_df = count_df.groupby(['activity', 'hash'], sort=True).size().reset_index(name='count')
    df = pd.merge(df, count_df, on=['activity', 'hash'], how='outer')
    df = df.sort_values(by=['activity', 'count', 'hash'], ascending=False)
    return df[['activity', 'hash', 'count', 'mol_name', 'isomer_id', 'conf_id', 'feature_ids']]


def strategy_extract_trainset(df, clust_strategy):
    if clust_strategy == 2:
        df = df.sort_values(by=['recall', 'F2', 'F05'], ascending=False).reset_index(drop=True)
        if df['F2'].iloc[0] == 1.0:
            df = df[(df['recall'] == 1.0) & (df['F2'] == 1.0)]
        elif df[df['F2'] >= 0.8].shape[0] <= 100:
            df = df[(df['recall'] == 1) & (df['F2'] >= 0.8)]
        else:
            df = df[(df['recall'] == 1) & (df['F2'] >= df['F2'].loc[100])]
    elif clust_strategy == 1:
        df = df.sort_values(by=['recall', 'F05', 'F2'], ascending=False).reset_index(drop=True)
        df = df[df['F05'] >= 0.8] if df[df['F05'] >= 0.8].shape[0] <= 100 else df[df['F05'] >= df['F05'].loc[100]]
    return df


# return type DataFrame: columns=['hash', 'TP', 'FP', 'precision', 'recall', 'F2', 'F05']
def calc_internal_stat(df, positives, clust_strategy, designating):
    if df[df['activity'] == designating[1]].empty:
        df['FP'] = [0] * df.shape[0]
        df = df.rename(columns={'count': 'TP'})
    else:
        df = df[df['activity'] == designating[0]].rename(columns={'count': 'TP'}).merge(
             df[df['activity'] == designating[1]].rename(columns={'count': 'FP'}),
             on='hash', how='right')
        df.loc[df['FP'].isnull(), 'FP'] = 0
    df = df[['hash', 'TP', 'FP']]
    df['precision'] = round(df['TP'] / (df['TP'] + df['FP']), 3)
    df['recall'] = round(df['TP'] / positives, 3)
    df['F2'] = round(5 * ((df['precision'] * df['recall']) / (4 * df['precision'] + df['recall'])), 3)
    df['F05'] = round(1.25 * ((df['precision'] * df['recall']) / (0.25 * df['precision'] + df['recall'])), 3)
    df = df[['hash', 'TP', 'FP', 'precision', 'recall', 'F2', 'F05']]
    # difference ways to check out the best models
    df = strategy_extract_trainset(df, clust_strategy)
    return df


def save_models_xyz(db, df_sub, path_pma, bin_step, cluster_id, num_ids):
    data = df_sub.drop_duplicates(subset=['hash']).values
    for num, (_, hash, count, mol_name, isomer_id, conf_id, feature_ids) in enumerate(data):
        pharm = Pharmacophore(bin_step=bin_step, cached=True)
        pharm.load_from_feature_coords(db.get_pharm(mol_name)[isomer_id][conf_id])
        pharm.save_to_xyz(os.path.join(path_pma, f'{cluster_id}_f{num_ids}_p{num}.xyz'),
                          tuple(map(int, feature_ids.split(','))))
    return len(data)


def gen_pharm_models(in_db, out_pma, trainset, tolerance, bin_step, current_nfeatures, upper, nfeatures, save_statistics):
    time_start = time.time()
    out_pma = os.path.join(out_pma, os.path.splitext(os.path.basename(in_db))[0])
    os.makedirs(out_pma, exist_ok=True)
    cluster_id = os.path.splitext(os.path.basename(trainset))[0]
    designating = ['1', '0']  # molecular activity
    clust_strategy = 1 if cluster_id == 'centroids' else 2
    positives = len([line for line in open(trainset).readlines() if line.strip().split()[2] == designating[0]])
    db = DB(in_db, flag='r')
    df_sub = gen_models(_gen_quadruplets(db, trainset, current_nfeatures, tolerance, bin_step))
    df = calc_internal_stat(df_sub[['activity', 'hash', 'count']].drop_duplicates(subset=['activity', 'hash']),
                            positives, clust_strategy, designating)
    if df.empty:
        sys.stderr.write(f'no {current_nfeatures}-points pharmacophore models for {cluster_id} training set\n')
        sys.exit(0)

    if save_statistics:
        path_files = os.path.join(out_pma, 'intermediate_statistics')
        os.makedirs(path_files, exist_ok=True)
        save_statistics = [path_files, cluster_id]

    df_sub = _keep_best_models(df, df_sub, save_statistics, current_nfeatures)
    if nfeatures is not None:
        if current_nfeatures >= nfeatures:
            _ = save_models_xyz(db, df_sub[df_sub['activity'] == designating[0]], out_pma, bin_step, cluster_id, current_nfeatures)

    while True:
        if current_nfeatures == upper:
            break
        current_nfeatures += 1
        df_sub_2 = gen_models(_plus_one_feature(db, df_sub, bin_step))
        df = calc_internal_stat(df_sub_2[['activity', 'hash', 'count']].drop_duplicates(subset=['activity', 'hash']),
                            positives, clust_strategy, designating)
        if df.empty:
            break

        df_sub = _keep_best_models(df, df_sub_2, save_statistics, current_nfeatures)
        if nfeatures is not None:
            if current_nfeatures >= nfeatures:
                _ = save_models_xyz(db, df_sub[df_sub['activity'] == designating[0]], out_pma, bin_step, cluster_id, current_nfeatures)

    num_models = save_models_xyz(db, df_sub[df_sub['activity'] == designating[0]], out_pma, bin_step, cluster_id, current_nfeatures)
    sys.stderr.write(f'train set {cluster_id}: {num_models} models ({round(time.time()-time_start, 3)}s)\n')
    sys.stderr.flush()


def create_parser():
    parser = argparse.ArgumentParser(description='Iteratively create ligand-based pharmacophore models.',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-db', '--database', metavar='active.db', required=True,
                        help='Input SQL database file with active compounds')
    parser.add_argument('-o', '--models', metavar='output/path', required=False, default=None,
                        help='Output path to a folder where will be saved the created pharmacophore models.'
                             'If omitted, the path will be generated automatically relative to project directory.')
    parser.add_argument('-ts', '--trainset', metavar='training_set.txt', required=True,
                        help='Path to tab-separated txt file with information about molecules from a training set.'
                             'Columns: SMILES, MOL_ID, ACTIVITY')
    parser.add_argument('-tol', '--tolerance', metavar='NUMERIC', type=float, default=0,
                        help='tolerance used for calculation of a stereoconfiguration sign.')
    parser.add_argument('-b', '--bin_step', type=int, required=False, default=1,
                        help='binning step.')
    parser.add_argument('-u', '--upper', type=int,  required=False, default=None,
                        help='limit the upper number of features in generated pharmacophores. '
                             'If omitted pharmacophores of maximum complexity will be generated.')
    parser.add_argument('-l', '--lower', type=int, default=3,
                        help='starting from this number of features, pharmacophore models will be created')
    parser.add_argument('-sm', '--save_model_complexity', type=int, required=False, default=None,
                        help='All pharmacophore models will be saved starting from this number of features.'
                             'If omitted will be saved only the most complex pharmacophore models')
    parser.add_argument('-s', '--save_statistics', type=bool, required=False, default=False,
                        help='If True the intermediate statistics of pharmacophore models generation will be stored')
    return parser


if __name__ == '__main__':
    parser = create_parser()
    args = parser.parse_args()
    gen_pharm_models(in_db=os.path.abspath(args.path_database),
                     trainset=os.path.abspath(args.trainset),
                     out_pma=os.path.abspath(args.path_models) if args.path_models else os.path.join(os.path.split(os.path.abspath(args.path_database))[0], 'models'),
                     bin_step=int(args.bin_step),
                     tolerance=args.tolerance,
                     current_nfeatures=args.lower,
                     upper=int(args.upper) if args.upper is not None else None,
                     nfeatures=int(args.feature) if args.upper is not None else None,
                     save_statistics=args.save_statistics)
