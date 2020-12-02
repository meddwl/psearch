#!/usr/bin/env python3
# author          : Alina Kutlushina
# date            : 02.12.20
# license         : BSD-3
# ==============================================================================

__author__ = 'Alina Kutlushina'

import os
import pandas as pd

import argparse


def calc_probability(df_vs, df_precision):
    df = df_vs.multiply(df_precision['precision'], axis=0)
    df.loc['probability'] = df.sum(axis=0)
    return df.loc['probability'].transpose()


def input_processing(pp_vs, list_mol_ids):
    if os.stat("file").st_size == 0:
        return 0
    list_pharms = [os.path.splitext(ff)[0] for ff in os.listdir(pp_vs)]
    df = pd.DataFrame(columns=list_mol_ids, index=list_pharms)
    for ff in os.listdir(pp_vs):
        ph = os.path.splitext(ff)[0]
        mols = [i.strip().split()[0] for i in open(os.path.join(pp_vs, ff)).readlines()]
        for mol_id in mols:
            df.at[ph, mol_id] = 1
    df.fillna(0, inplace=True)
    return df


def main(pp_vs, pp_smi, pp_models_stat, pp_output):
    mol_ids = [i.strip().split()[1] for i in open(pp_smi).readlines()]
    df_models_stat = pd.read_csv(pp_models_stat, sep='\t', index_col='model')
    df_vs = input_processing(pp_vs, mol_ids)
    df_res = calc_probability(df_vs, df_models_stat)
    df_res.to_csv(pp_output, sep='\t')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Determination of the probability of activity of a molecule(-s)'
                                                 'based on pharmacophore VS result',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-vs', '--path_vs', metavar='path/to/vs/res', required=True,
                        help='path to the virtual screening result')
    parser.add_argument('-m', '--molecules', metavar='molecules.smi', required=True,
                        help='.smi file containing the molecules smiles and their ids')
    parser.add_argument('-p', '--models_stat', metavar='pharmacophores_stat.csv', required=True,
                        help='.csv file with the precision of pharmacophore models')
    parser.add_argument('-o', '--output', metavar='external_statistics.txt', default=None,
                        help='output text file where will be saved the prediction')

    args = parser.parse_args()
    main(args.path_vs, args.molecules, args.models_stat, args.output)
