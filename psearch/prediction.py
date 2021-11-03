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
default_modelstat = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'pharmacophores', 'pharmacophores_stat.csv')

def calc_probability(df_vs, df_precision, target_id, scoring_scheme):
    df = df_vs.mul(df_precision['precision'].astype(float), axis=0)
    if scoring_scheme == 'mean':
        df.loc[target_id] = df.mean(axis=0, skipna=True)
    elif scoring_scheme == 'max':
        df.loc[target_id] = df.max(axis=0, skipna=True)
    return round(df.loc[target_id], 3)


def input_processing(list_vs, models_list):
    df = pd.DataFrame(index=models_list)
    for ff in list_vs:
        ph = os.path.splitext(os.path.basename(ff))[0].split('.')[1]
        mols = [i.strip().split()[0] for i in open(ff).readlines()]
        for mol_id in mols:
            df.at[ph, mol_id] = 1
    return df


def main(pp_vs, pp_models_stat, scoring_scheme, pp_output):
    df_models_stat = pd.read_csv(pp_models_stat, sep='\t', index_col='model_id')
    vs_files = defaultdict(list)
    for fname in os.listdir(pp_vs):
        target_id = os.path.splitext(fname)[0].split('.')[0]
        vs_files[target_id].append(os.path.join(pp_vs, fname))

    res = pd.DataFrame()
    for target_id, fvs in vs_files.items():
        df_models = df_models_stat[df_models_stat['target_id'] == target_id]
        df_vs = input_processing(fvs, df_models.index.tolist())
        df_res = calc_probability(df_vs, df_models, target_id, scoring_scheme)
        res = res.merge(df_res, left_index=True, right_index=True, how='outer')
    # res = res.transpose()
    res.index.name = 'mol_id'
    res.sort_values(by=res.columns.tolist()[0], ascending=False, inplace=True)
    res.to_csv(pp_output, sep='\t')


def entry_point():
    parser = argparse.ArgumentParser(description='Determination of the probability of activity of a molecule(-s)'
                                                 'based on pharmacophore VS result',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-vs', '--path_vs', metavar='path/to/vs/res', required=True,
                        help='path to the virtual screening result')
    parser.add_argument('-stat', '--models_stat', metavar='pharmacophores_stat.csv',
                        default=default_modelstat,
                        help='file with the calculated precision of pharmacophore models. '
                             'By default, statistics are used for the chembl models, which are in pserch.')
    parser.add_argument('-s', '--scoring_scheme', default='mean',
                        help='two schemes (Max and Mean) of probability calculation for consensus prediction '
                             'based on individual pharmacophore models were proposed')
    parser.add_argument('-o', '--output', metavar='external_statistics.txt', default=None,
                        help='output text file where will be saved the prediction')

    args = parser.parse_args()
    main(os.path.abspath(args.path_vs), os.path.abspath(args.models_stat),
         args.scoring_scheme, os.path.abspath(args.output))


if __name__ == '__main__':
    entry_point()
