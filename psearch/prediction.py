#!/usr/bin/env python3
# author          : Alina Kutlushina
# date            : 02.12.20
# license         : BSD-3
# ==============================================================================

__author__ = 'Alina Kutlushina'

import os
import pandas as pd
import argparse
from collections import defaultdict
default_pharmstat = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'pharmacophores', 'pharmacophores_stat.csv')


def calc_probability(df_vs, df_precision, target_id, scoring_scheme):
    df = df_vs.mul(df_precision['precision'].astype(float), axis=0)
    if scoring_scheme == 'mean':
        df.loc[target_id] = df.mean(axis=0, skipna=True)
    elif scoring_scheme == 'max':
        df.loc[target_id] = df.max(axis=0, skipna=True)
    return round(df.loc[target_id], 3)


def input_processing(list_vs, pharm_list):
    df = pd.DataFrame(index=pharm_list)
    for ff in list_vs:
        pharm_id = os.path.splitext(os.path.basename(ff))[0].split('.')[1]
        mols_id = [i.strip().split()[0] for i in open(ff).readlines()]
        df.loc[pharm_id, mols_id] = 1
    return df


def main(pp_vs, pp_pharm_stat, scoring_scheme, min_num_features, pp_output):
    df_pharm_stat = pd.read_csv(pp_pharm_stat, sep='\t', index_col='model_id')
    vs_files = defaultdict(list)

    for fname in os.listdir(pp_vs):
        num_features = int(fname.split('.')[1].split('_')[1].split('pharm')[1])  # we need to find another solution how to check pharmacophore complexity
        if num_features >= min_num_features:
            target_id = os.path.splitext(fname)[0].split('.')[0]
            vs_files[target_id].append(os.path.join(pp_vs, fname))

    res = pd.DataFrame()
    target_ids = sorted(vs_files.keys())

    for target_id in target_ids:
        df_pharm = df_pharm_stat[df_pharm_stat['target_id'] == target_id]
        df_vs = input_processing(vs_files[target_id], df_pharm.index.tolist())
        df_res = calc_probability(df_vs, df_pharm, target_id, scoring_scheme)
        res = res.merge(df_res, left_index=True, right_index=True, how='outer')

    res.index.name = 'cmp'
    res.sort_values(by=res.columns.tolist()[0], ascending=False, inplace=True)
    res.to_csv(pp_output, sep='\t')


def entry_point():
    parser = argparse.ArgumentParser(description='Determination of the probability of activity of a molecule(-s)'
                                                 'based on pharmacophore VS result',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--path_vs', metavar='DIRNAME', required=True,
                        help='path to the virtual screening result')
    parser.add_argument('-p', '--pharm_stat', metavar='FILENAME', default=None,
                        help='file with the calculated precision of pharmacophores. '
                             'By default, statistics of psearch pharmacophores are used.'
                             'Required headers: "target_id", "model_id", "precision"')
    parser.add_argument('-s', '--scoring_scheme', metavar='KEYWORD', default='mean',
                        help='two schemes (Max and Mean) of probability calculation for consensus prediction '
                             'based on individual pharmacophores were proposed')
    parser.add_argument('-f', '--min_features', metavar='INTEGER', default=None, type=int,
                            help='minimum number of features with distinct coordinates in models. Models having less '
                                 'number of features will be skipped. Default: all models will be screened.')
    parser.add_argument('-o', '--output', metavar='FILENAME', default=None,
                        help='output text file where will be saved the prediction')

    args = parser.parse_args()
    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    pharm_stat = os.path.abspath(args.pharm_stat) if args.pharm_stat else default_pharmstat

    main(pp_vs=os.path.abspath(args.path_vs),
         pp_pharm_stat=pharm_stat,
         scoring_scheme=args.scoring_scheme,
         min_num_features=args.min_features if args.min_features != None else 0,
         pp_output=os.path.abspath(args.output))

if __name__ == '__main__':
    entry_point()
