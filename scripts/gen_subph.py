#!/usr/bin/env python3
# author          : Alina Kutlushina
# date            : 01.05.2018
# license         : BSD-3
#==============================================================================

import os
import sys
import time
import argparse
import numpy as np
import pandas as pd
import sqlite3 as sql
from collections import defaultdict
from pmapper.pharmacophore import Pharmacophore


def get_pharm(cur, conf_id):
    cur.execute('SELECT feature_label, x, y, z FROM feature_coords WHERE conf_id = ?', (conf_id,))
    res = cur.fetchall()
    if res:
        p = Pharmacophore()
        p.load_from_feature_coords(tuple((i[0], tuple(i[1:])) for i in res))
        return p
    else:
        return None


def extract_better_models(df, cur, out, nids):
    start_time = time.time()
    d = defaultdict(list)
    for c, l, m in zip(df['conf_id'], df['feature_ids'], df['mol_name']):
        d[c, m].append(tuple(int(i) for i in l.split(',')))

    d_nids = defaultdict(list)
    for (conf_id, mol_name), ids in d.items():
        p = get_pharm(cur, int(conf_id))
        if p:
            for hash, labels in p.iterate_pharm1(ids):
                d_nids['hash'].append(hash)
                d_nids['mol_name'].append(mol_name)
                d_nids['conf_id'].append(conf_id)
                d_nids['feature_ids'].append(','.join('{}'.format(i) for i in labels))

    df_nids = pd.DataFrame(d_nids)
    df_nids = df_nids.drop_duplicates(subset=['mol_name', 'hash'])
    count_df = df_nids.groupby(['hash'], sort=True).size().reset_index(name='count')
    count_df = count_df.sort_values(by=['count'], ascending=False)
    df = pd.merge(count_df, df_nids, on='hash', how='left')
    df = df.sort_values(by=['count', 'hash'], ascending=False)
    df = df[['hash', 'count', 'mol_name', 'conf_id', 'feature_ids']]
    df.to_csv(out, index=None, sep='\t')
    sys.stderr.write('{}: ({}s)\n'.format(os.path.basename(out), round(time.time() - start_time, 3)))
    if df.empty:
        return None
    else:
        return df


def strategy_extract_trainset(df, make_clust):
    df = df.reset_index(drop=True) # re-indexing
    if make_clust:
        df = df.sort_values(by=['recall', 'F2', 'F05'], ascending=False)
        if df['F2'].iloc[0] == 1.0:
            df = df[(df['recall'] == 1.0) & (df['F2'] == 1.0)]
        elif df[df['F2'] >= 0.8].shape[0] <= 100:
            df = df[(df['recall'] == 1) & (df['F2'] >= 0.8)]
        else:
            df = df[(df['recall'] == 1) & (df['F2'] >= df['F2'].loc[100])]

    else:
        df = df.sort_values(by=['recall', 'F05', 'F2'], ascending=False)
        if df[df['F05'] >= 0.8].shape[0] <= 100:
            df = df[df['F05'] >= 0.8]
        else:
            df = df[df['F05'] >= df['F05'].loc[100]]
    return df


def save_int_stat(df_act, df_inact, act_trainset, out, make_clust, cluster_num, num):
    #
    start_time = time.time()
    df_act = df_act.rename(columns={'count': 'TP'})
    df_act = df_act.drop(['mol_name', 'conf_id', 'feature_ids'], axis=1)
    df_act = df_act.drop_duplicates(subset=['hash'])
    
    if not df_inact.empty:
        df_inact = df_inact.rename(columns={'count': 'FP'})
        df_inact = df_inact.drop(['mol_name', 'conf_id', 'feature_ids'], axis=1)
        df_inact = df_inact.drop_duplicates(subset=['hash'])
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
    df['model'] = ['{}_pharm{}'.format(cluster_num, num) + '_%i' % x for x in range(df.shape[0])]
    df = df[['model', 'hash', 'TP', 'FP', 'precision', 'recall', 'F2', 'F05']]
    # difference ways to check out the best models
    df = strategy_extract_trainset(df, make_clust)
    if not df.empty:
        df.to_csv(out, index=None, sep='\t')
    sys.stderr.write('{}: ({}s)\n'.format(os.path.basename(out), round(time.time() - start_time, 3)))
    return df


def exrtact_subph_pma(df, cur, out_pma):
    time_start = time.time()
    df = df.drop_duplicates(subset=['hash'])
    i = 0
    for conf_id, labels, name_model in zip(df['conf_id'], df['feature_ids'], df['model']):
        p = get_pharm(cur, int(conf_id))
        if p:
            i += 1
            p.save_to_pma(os.path.join(out_pma, name_model + '.pma'), tuple(map(int, labels.split(','))))
    sys.stderr.write('{}: extract {} models passed ({}s)\n\n'.format(os.path.split(out_pma)[1],
                                                                     i, round(time.time() - time_start, 3)))


def main(sub_act, sub_inact, in_adb, in_indb, act_trainset, lower):
    """
    :param in_internstat: input txt file with internal statistics
    :param sub_act: input txt file with active subpharmacophore models
    :param sub_inact: input txt file with inactive subpharmacophore models
    :param in_adb: input SQL databased with active compounds
    :param in_indb: input SQL databased with inactive compounds
    :param lower: minimum (starting) number of nodes of the subpharmacophore models
    :return: pma files best models and txt files with internal statistics
    """

    def make_dir(new_path):
        if not os.path.exists(new_path):
            os.mkdir(new_path)

    def path_out_subph(subph, belonging, tr_num, num):
        return os.path.join(os.path.dirname(subph), 'ph_{}_{}_pharm{}.txt'.format(belonging, tr_num, num))

    def path_out_interstat(in_adb, tr_num, num):
        return os.path.join(os.path.dirname(in_adb), 'internal_statistics_{}_pharm{}.txt'.format(tr_num, num))

    if os.path.basename(act_trainset).split('.')[0].split('_')[1] == 'centroid':
        make_clust = False
    else:
        make_clust = True

    tr_num = os.path.basename(sub_act).split('.')[0].split('_')[2]
    num = lower

    df_act = pd.read_csv(sub_act, sep='\t')
    df_inact = pd.read_csv(sub_inact, sep='\t')
    
    out_inetsts = path_out_interstat(in_adb, tr_num, num)
    df = save_int_stat(df_act, df_inact, act_trainset, out_inetsts, make_clust, tr_num, num)

    if df.empty:
        return None, 0

    df = df.drop(['TP', 'FP', 'precision', 'recall', 'F2', 'F05'], axis=1)
    df_act = pd.merge(df, df_act, on='hash', how='left')
    df_inact = pd.merge(df, df_inact, on='hash', how='inner')
    
    with sql.connect(in_adb) as conn_adb:
        cur_adb = conn_adb.cursor()

        out_pma = os.path.join(os.path.dirname(os.path.abspath(in_adb)), 'models', tr_num, 'pharm%i' % lower)
        make_dir(os.path.split(os.path.split(out_pma)[0])[0])
        make_dir(os.path.split(out_pma)[0])
        make_dir(out_pma)

        exrtact_subph_pma(df_act, cur_adb, out_pma)

        while not df.empty:

            num += 1
            print(num)
            out_ph_act = path_out_subph(sub_act, 'active', tr_num, num)
            df_act = extract_better_models(df_act, cur_adb, out_ph_act, num)

            out_inetsts = path_out_interstat(in_adb, tr_num, num)

            if not df_inact.empty:
                with sql.connect(in_indb) as conn_indb:
                    cur_indb = conn_indb.cursor()
                    out_ph_inact = path_out_subph(sub_inact, 'inactive', tr_num, num)
                    df_inact = extract_better_models(df_inact, cur_indb, out_ph_inact, num)
                    df = save_int_stat(df_act, df_inact, act_trainset, out_inetsts, make_clust, tr_num, num)
            else:
                df = save_int_stat(df_act, df_inact, act_trainset, out_inetsts, make_clust, tr_num, num)

            if not df.empty:
                out_pma = os.path.join(os.path.split(out_pma)[0], 'pharm%i' % num)
                make_dir(out_pma)
                df_act = pd.merge(df, df_act.drop(['count'], axis=1), on=['hash'], how='left')
                exrtact_subph_pma(df_act, cur_adb, out_pma)
            else:
                sys.stderr.write('{}\n'.format(os.path.basename(out_inetsts)))
                return out_pma, num-1
                        

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='ilter out duplicated hash-stereo pair for each compound ID (without stereo)')
    parser.add_argument('-sub_act', '--in_subph_active', metavar='ph_active_trYY_pharmXX.txt', required=True,
                        help='')
    parser.add_argument('-sub_inact', '--in_subph_inactive', metavar='ph_active_trYY_pharmXX.txt', required=True,
                        help='')
    parser.add_argument('-in_adb', '--in_active_database', metavar='active.db', required=True,
                        help='input SQL database file with active compounds')
    parser.add_argument('-in_indb', '--in_inactive_database', metavar='inactive.db', required=True,
                        help='input SQL database file with active compounds')
    parser.add_argument('-iats', '--in_active_trainset', metavar='ph_inactives_subsetY_pX.txt', required=True,
                        help='txt file with information adout active models: '
                             'model, hash, stereo, nact, ninact, nact/ninact, conf_id, feature_ids')
    parser.add_argument('-clust', '--make_clust', action='store_true', default=False,
                        help='if set training sets will be created for separate clusters, '
                             'otherwise only one training set will be created.')
    parser.add_argument('-nids', '--number_ids', default=0,
                        help='number of features of input models')

    args = vars(parser.parse_args())
    for o, v in args.items():
        if o == "in_subph_active": sub_act = v
        if o == "in_subph_inactive": sub_inact = v
        if o == "in_active_database": in_adb = v
        if o == "in_inactive_database": in_indb = v
        if o == "in_active_trainset": act_trainset = v
        if o == "number_ids": nids = int(v)

    main(sub_act=sub_act,
         sub_inact=sub_inact,
         in_adb=in_adb,
         in_indb=in_indb,
         act_trainset=act_trainset,
         lower=nids)
