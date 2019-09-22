#!/usr/bin/env python3
# author          : Alina Kutlushina
# date            : 01.05.2018
# license         : BSD-3
#==============================================================================

import os, sys
import time
import json
import pandas as pd
import argparse
from pmapper.pharmacophore import Pharmacophore as P


def max_edge(in_model):
    model = os.path.abspath(in_model)
    p = P()
    p.load_from_pma(model)
    coords = p.get_feature_coords()

    edge = 0
    for i, c1 in enumerate(coords):
        for j, c2 in enumerate(coords[i + 1:]):
            e = ((c1[1][0] - c2[1][0]) ** 2 + (c1[1][1] - c2[1][1]) ** 2 + (c1[1][2] - c2[1][2]) ** 2) ** (1/2)
            if e > edge:
                edge = e
    return edge


def get_external_stat(mol_act, mol_inact, in_pma, in_act_screen, in_inact_screen):

    medge = max_edge(in_pma)
    model = os.path.basename(in_pma).split('.')[0]
    with open(in_pma) as fpma:
        d = json.loads(fpma.readline().strip())
        labels = ''.join(i[0] for i in d['feature_coords'])
        num_uniq_features = len(set(tuple(feature[1]) for feature in d['feature_coords']))
        
    all_active_mol = len(open(mol_act).readlines())
    all_inactive_mol = len(open(mol_inact).readlines())
    if not os.path.exists(in_act_screen):
        TP = 0
    else:
        TP = len(open(in_act_screen).readlines())
    if not os.path.exists(in_inact_screen):
        FP = 0
    else:
        FP = len(open(in_inact_screen).readlines())
    FN = all_active_mol - TP
    TN = all_inactive_mol - FP

    try:
        precision = TP / (TP + FP)
        recall = TP / (TP + FN)
        FPR = FP / (TN + FP)
        spc = TN / (TN + FP)
        f1 = (2 * precision * recall) / (precision + recall)
        f2 = (5 * precision * recall) / (4 * precision + recall)
        f05 = (1.25 * precision * recall) / (0.25 * precision + recall)
        ba = round((recall + spc) / 2, 3)
        ef = (TP / (TP + FP)) / (all_active_mol / (all_inactive_mol + all_active_mol))

        return [model,
                TP,
                FP,
                round(precision, 3),
                round(FPR, 3),
                round(recall, 3),
                round(f1, 3),
                round(f2, 3),
                round(f05, 3),
                round(ba, 3),
                round(ef, 3),
                num_uniq_features,
                medge,
                labels]
    
    except ZeroDivisionError:
        return [model,
                TP,
                FP,
                '-',
                '-',
                '-',
                '-',
                '-',
                '-',
                '-',
                '-',
                num_uniq_features,
                medge,
                labels]


def main(mol_act, mol_inact, path_to_pma, in_act_screen, in_inact_screen, out_external):
    
    start_time = time.time()

    df_result = pd.DataFrame(columns=['model', 'TP', 'FP', 'precision', 'FPR', 'recall',
                                      'F1', 'F2', 'F05', 'BA', 'EF', 'num_uniq_F', 'max_edge', 'features'])


    for enum, in_pma in enumerate(os.listdir(path_to_pma)):
        ppath = (os.path.join(os.path.abspath(in_act_screen), '{}.txt'.format(in_pma.split('.')[0])),
                 os.path.join(os.path.abspath(in_inact_screen), '{}.txt'.format(in_pma.split('.')[0])),
                 os.path.abspath(os.path.join(path_to_pma, in_pma)))

        result = get_external_stat(mol_act, mol_inact, ppath[2], ppath[0], ppath[1])
        if result:
            df_result.loc[enum] = result
        else:
            continue

    df_result = df_result.sort_values(by=['recall', 'F05', 'F2'], ascending=False)

    df_result.to_csv(out_external, index=None, sep='\t')
    sys.stderr.write('{}: ({}s)\n\n'.format(os.path.basename(out_external), round(time.time()-start_time,3)))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-ma', '--input_active_mol', metavar='active.smi', required=True,
                        help='.smi file with active molecules.')
    parser.add_argument('-mi', '--input_inactive_mol', metavar='inactive.smi', required=True,
                        help='.smi file with inactive molecules.')
    parser.add_argument('-p', '--path_to_pma', metavar='models/pYY/', required=True,
                        help='path to pma files')
    parser.add_argument('-as', '--act_screen', metavar='screen/pYY/', required=True,
                        help='path to screen act')
    parser.add_argument('-is', '--inact_screen', metavar='screen/pYY/',
                        help='path to screen inact')
    parser.add_argument('-o', '--out_external', metavar='external_subsetXX_pYY.txt', default=None,
                        help='output text file, which will contain: model, n_act, n_inact, precision, recall, '
                             'F1, BA, EF, features, n_act_train, n_inact_train. ')

    args = vars(parser.parse_args())
    for o, v in args.items():
        if o == "input_active_mol": mol_act = v
        if o == "input_inactive_mol": mol_inact = v
        if o == "path_to_pma": path_to_pma = v
        if o == "act_screen": act_screen = v
        if o == "inact_screen": inact_screen = v
        if o == "out_external": out_external = v

    if out_external is None:
        out_external = os.path.join(os.path.split(os.path.dirname(mol_act))[0], 'result.txt')

    main(mol_act, mol_inact, path_to_pma, act_screen, inact_screen, out_external)
    
