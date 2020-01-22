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


def create_parser():
    parser = argparse.ArgumentParser(description='', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-ma', '--input_active', metavar='active.smi', required=True,
                        help='.smi file with active molecules.')
    parser.add_argument('-mi', '--input_inactive', metavar='inactive.smi', required=True,
                        help='.smi file with inactive molecules.')
    parser.add_argument('-t', '--path_trainset', metavar='active.smi', required=True,
                        help='path to folders with files of training sets.')
    parser.add_argument('-p', '--path_to_pma', metavar='models/pYY/', required=True,
                        help='path to pma files')
    parser.add_argument('-as', '--act_screen', metavar='screen/pYY/', required=True,
                        help='path to screen act')
    parser.add_argument('-is', '--inact_screen', metavar='screen/pYY/', required=True,
                        help='path to screen inact')
    parser.add_argument('-o', '--out_external', metavar='external_subsetXX_pYY.txt', default=None,
                        help='output text file, which will contain next columns: model, n_act, n_inact, precision, '
                             'recall, F1, BA, EF, features, n_act_train, n_inact_train. ')
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


def get_external_stat(mol_act, mol_inact, ts_act, ts_inact, in_pma, in_act_screen, in_inact_screen):
    medge = max_edge(in_pma)
    model = os.path.basename(in_pma).split('.')[0]
    with open(in_pma) as fpma:
        d = json.loads(fpma.readline().strip())
        labels = ''.join(i[0] for i in d['feature_coords'])
        num_uniq_features = len(set(tuple(feature[1]) for feature in d['feature_coords']))

    ts_active_mol = [ii.strip().split('\t')[1] for ii in open(ts_act).readlines()]
    if os.path.exists(in_act_screen):
        act_screen = [ii.strip().split('\t')[0] for ii in open(in_act_screen).readlines()]
        act_screen = set(act_screen).difference(ts_active_mol)
    else:
        act_screen = []
    ts_inactive_mol = [ii.strip().split('\t')[1] for ii in open(ts_inact).readlines()]
    if os.path.exists(in_inact_screen):
        inact_screen = [ii.strip().split('\t')[0] for ii in open(in_inact_screen).readlines()]
        inact_screen = set(inact_screen).difference(ts_inactive_mol)
    else:
        inact_screen = []

    p = len(open(mol_act).readlines()) - len(ts_active_mol)
    n = len(open(mol_inact).readlines()) - len(ts_inactive_mol)
    tp = len(act_screen)
    fp = len(inact_screen)
    fn = p - tp
    tn = n - fp

    if tp == 0 and fp == 0:
        return [model, tp, fp, p, n,
                np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan,
                num_uniq_features, medge, labels]

    elif tp != 0 or fp != 0:
        precision = round(tp / (tp + fp), 3)
        recall = round(tp / (tp + fn), 3)
        fpr = round(fp / (tn + fp), 3)
        ba = round((recall + (tn / (tn + fp))) / 2, 3)
        ef = round((tp / (tp + fp)) / (p / (p + n)), 3)
        if precision != 0 and recall != 0:
            f1 = round((2 * precision * recall) / (precision + recall), 3)
            f2 = round((5 * precision * recall) / (4 * precision + recall), 3)
            f05 = round((1.25 * precision * recall) / (0.25 * precision + recall), 3)
        else:
            f1 = np.nan
            f2 = np.nan
            f05 = np.nan
        return [model, tp, fp, p, n, precision, fpr, recall, f1, f2, f05, ba, ef, num_uniq_features, medge, labels]


def calc_stat(mol_act, mol_inact, path_ts, path_to_pma, in_act_screen, in_inact_screen, out_external):
    start_time = time.time()

    df_result = pd.DataFrame(columns=['model', 'TP', 'FP', 'P', 'N', 'precision', 'FPR', 'recall',
                                      'F1', 'F2', 'F05', 'BA', 'EF', 'num_uniq_F', 'max_edge', 'features'])

    for enum, in_pma in enumerate(os.listdir(path_to_pma)):
        ppath = (os.path.join(os.path.abspath(path_ts), 'active_{}.csv'.format(in_pma.split('_')[0])),
                 os.path.join(os.path.abspath(path_ts), 'inactive_{}.csv'.format(in_pma.split('_')[0])),
                 os.path.join(os.path.abspath(in_act_screen), '{}.txt'.format(in_pma.split('.')[0])),
                 os.path.join(os.path.abspath(in_inact_screen), '{}.txt'.format(in_pma.split('.')[0])),
                 os.path.abspath(os.path.join(path_to_pma, in_pma)))

        result = get_external_stat(mol_act, mol_inact, ppath[0], ppath[1], ppath[4], ppath[2], ppath[3])
        if result:
            df_result.loc[enum] = result
        else:
            continue

    df_result = df_result.sort_values(by=['recall', 'F05', 'F2'], ascending=False)

    df_result.to_csv(out_external, index=None, sep='\t')
    sys.stderr.write('{}: ({}s)\n\n'.format(os.path.basename(out_external), round(time.time() - start_time, 3)))


def entry_point():
    parser = create_parser()
    args = vars(parser.parse_args())
    for o, v in args.items():
        if o == "input_active": mol_act = v
        if o == "input_inactive": mol_inact = v
        if o == "path_trainset": path_ts = v
        if o == "path_to_pma": path_to_pma = v
        if o == "act_screen": act_screen = v
        if o == "inact_screen": inact_screen = v
        if o == "out_external": out_external = v

    if out_external is None:
        out_external = os.path.join(os.path.split(os.path.dirname(mol_act))[0], 'result.txt')

    if not os.path.exists(os.path.dirname(out_external)):
        os.makedirs(os.path.dirname(out_external))

    calc_stat(mol_act, mol_inact, path_ts, path_to_pma, act_screen, inact_screen, out_external)


if __name__ == '__main__':
    entry_point()
