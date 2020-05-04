#!/usr/bin/env python3
# author          : Alina Kutlushina
# date            : 01.05.2018
# license         : BSD-3
# ==============================================================================

import os
import sys
import time
import json
import pandas as pd
import numpy as np
import argparse
from pmapper.pharmacophore import Pharmacophore as P
from rdkit import Chem


def create_parser():
    parser = argparse.ArgumentParser(description='calculate external statistics',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-m', '--path_molecules', metavar='molecules.smi', required=True,
                        help='.smi file with active and inactive molecules.')
    parser.add_argument('-t', '--path_trainset', metavar='path/to/trainset', required=True,
                        help='path to folders with files of training sets.')
    parser.add_argument('-p', '--path_to_pma', metavar='path/to/models', required=True,
                        help='path to pma files')
    parser.add_argument('-s', '--path_screen', metavar='path/to/screen', required=True,
                        help='path to screen results')
    parser.add_argument('-o', '--out_external', metavar='external_statistics.txt', default=None,
                        help='output text file where will be saved external statistics')
    return parser


def max_edge(in_model):
    model = os.path.abspath(in_model)
    p = P()
    p.load_from_pma(model)
    coords = p.get_feature_coords()
    edge = 0
    for i, c1 in enumerate(coords):
        for j, c2 in enumerate(coords[i + 1:]):
            e = ((c1[1][0] - c2[1][0]) ** 2 + (c1[1][1] - c2[1][1]) ** 2 + (c1[1][2] - c2[1][2]) ** 2) ** (1 / 2)
            if e > edge:
                edge = e
    return edge


def get_external_stat(path_mols, ts_act, ts_inact, in_pma, in_screen):
    medge = max_edge(in_pma)
    model = os.path.splitext(os.path.basename(in_pma))[0]
    with open(in_pma) as fpma:
        d = json.loads(fpma.readline().strip())
        labels = ''.join(i[0] for i in d['feature_coords'])
        num_uniq_features = len(set(tuple(feature[1]) for feature in d['feature_coords']))

    ts_act_mol = [ii.strip().split()[1] for ii in open(ts_act).readlines()]
    ts_inact_mol = [ii.strip().split()[1] for ii in open(ts_inact).readlines()]

    df_mols = pd.read_csv(path_mols, sep='\t', header=None)
    df_mols.rename(columns={0: 'smiles', 1: 'mol_name', 2: 'activity'}, inplace=True)
    if not Chem.MolFromSmiles(df_mols.at[0, 'smiles']):
        df_mols.drop(index=0)
    if df_mols['activity'].dtypes == 'int64':
        df_act = df_mols[(df_mols['activity'] == 1) & (~df_mols['mol_name'].isin(ts_act_mol))]
        df_inact = df_mols[(df_mols['activity'] == 0) & (~df_mols['mol_name'].isin(ts_inact_mol))]
    else:
        df_act = df_mols[(df_mols['activity'] == 'active') & (~df_mols['mol_name'].isin(ts_act_mol))]
        df_inact = df_mols[(df_mols['activity'] == 'inactive') & (~df_mols['mol_name'].isin(ts_inact_mol))]

    if os.path.exists(in_screen):
        res_screen = [ii.strip().split()[0] for ii in open(in_screen).readlines()]
        act_screen = set(res_screen) & set(df_act['mol_name'])
        inact_screen = set(res_screen) & set(df_inact['mol_name'])
    else:
        act_screen = []
        inact_screen = []

    p = df_act.shape[0]
    n = df_inact.shape[0]
    tp = len(act_screen)
    fp = len(inact_screen)
    # fn = p - tp
    tn = n - fp

    recall = tp / p
    tnr = tn / n
    fpr = fp / n
    ba = (recall + tnr) / 2
    # accuracy = (tp + tn) / (p + n)

    if tp == 0 and fp == 0:
        precision = -1
        ef = -1
        f1 = -1
        f2 = -1
        f05 = -1
    else:
        precision = tp / (tp + fp)
        ef = precision / (p / (p + n))
        f1 = (2 * precision * recall) / (precision + recall)
        f2 = (5 * precision * recall) / (4 * precision + recall)
        f05 = (1.25 * precision * recall) / (0.25 * precision + recall)

    return [model, tp, fp, p, n, precision, fpr, recall, f1, f2, f05, ba, ef, num_uniq_features, medge, labels]


def calc_stat(path_mols, path_ts, path_pma, path_screen, out_external):
    start_time = time.time()
    if not os.path.exists(os.path.dirname(out_external)):
        os.makedirs(os.path.dirname(out_external))

    df_result = pd.DataFrame(columns=['model', 'TP', 'FP', 'P', 'N', 'precision', 'FPR', 'recall',
                                      'F1', 'F2', 'F05', 'BA', 'EF', 'num_uniq_F', 'max_edge', 'features'])

    for enum, in_pma in enumerate(os.listdir(path_pma)):
        ppath = (os.path.join(path_ts, f'active_{in_pma.split("_")[0]}.smi'),
                 os.path.join(path_ts, f'inactive_{in_pma.split("_")[0]}.smi'),
                 os.path.join(path_screen, f'{os.path.splitext(in_pma)[0]}.txt'),
                 os.path.abspath(os.path.join(path_pma, in_pma)))
        result = get_external_stat(path_mols, ppath[0], ppath[1], ppath[3], ppath[2])
        if result:
            df_result.loc[enum] = result
        else:
            continue

    df_result = df_result.sort_values(by=['recall', 'F05', 'F2'], ascending=False)
    df_result = round(df_result, 3)
    df_result.to_csv(out_external, index=None, sep='\t')
    sys.stderr.write(f'{os.path.basename(out_external)}: ({round(time.time() - start_time, 3)}s)\n\n')


if __name__ == '__main__':
    parser = create_parser()
    args = vars(parser.parse_args())
    for o, v in args.items():
        if o == "path_molecules": path_mols = os.path.abspath(v)
        if o == "path_trainset": path_ts = os.path.abspath(v)
        if o == "path_to_pma": path_to_pma = os.path.abspath(v)
        if o == "path_screen": path_screen = os.path.abspath(v)
        if o == "out_external": out_external = os.path.abspath(v)

    if out_external is None:
        out_external = os.path.join(os.path.dirname(path_mols), 'result.txt')

    calc_stat(path_mols, path_ts, path_to_pma, path_screen, out_external)
