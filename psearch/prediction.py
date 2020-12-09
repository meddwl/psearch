#!/usr/bin/env python3
# author          : Alina Kutlushina
# date            : 02.12.20
# license         : BSD-3
# ==============================================================================

__author__ = 'Alina Kutlushina'

import os
import pandas as pd
import argparse


def calc_probability(df_vs, df_precision, target_id):
    df = df_vs.multiply(df_precision['precision'], axis=0)
    df.loc[target_id] = df.sum(axis=0)
    return round(df.loc[target_id].transpose(), 3)


def input_processing(pp_vs, target_id, models_list):
    pp_vs_t = os.path.join(pp_vs, target_id)
    df = pd.DataFrame(index=models_list)
    if len(os.listdir(pp_vs_t)) == 0:
        return df
    for ff in os.listdir(pp_vs_t):
        ph = os.path.splitext(ff)[0]
        mols = [i.strip().split()[0] for i in open(os.path.join(pp_vs_t, ff)).readlines()]
        for mol_id in mols:
            df.at[ph, mol_id] = 1
    df.fillna(0, inplace=True)
    return df


def main(pp_vs, pp_models_stat, pp_output):
    df_models_stat = pd.read_csv(pp_models_stat, sep='\t', index_col='model')
    res = pd.DataFrame()
    for target_id in os.listdir(pp_vs):
        df_models = df_models_stat[df_models_stat['protein'] == target_id]
        df_vs = input_processing(pp_vs, target_id, df_models.index.tolist())
        if df_vs.empty:
            continue
        df_res = calc_probability(df_vs, df_models, target_id)
        res = res.merge(df_res, left_index=True, right_index=True, how='outer')
    res.index.name = 'pharmacophores_id'
    res.to_csv(pp_output, sep='\t')


def entry_point():
    parser = argparse.ArgumentParser(description='Determination of the probability of activity of a molecule(-s)'
                                                 'based on pharmacophore VS result',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-vs', '--path_vs', metavar='path/to/vs/res', required=True,
                        help='path to the virtual screening result')
    parser.add_argument('-p', '--models_stat', metavar='pharmacophores_stat.csv', required=True,
                        help='.csv file with the precision of pharmacophore models')
    parser.add_argument('-o', '--output', metavar='external_statistics.txt', default=None,
                        help='output text file where will be saved the prediction')

    args = parser.parse_args()
    main(args.path_vs, args.models_stat, args.output)


if __name__ == '__main__':
    entry_point()
