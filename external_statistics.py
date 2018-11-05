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


def get_external_stat(mol_act, mol_inact, ts_act, ts_inact, in_pma, in_act_screen, in_inact_screen):

    act_screen = []
    inact_screen = []

    ts_active_mol = []
    ts_inactive_mol = []
    all_active_mol = len(open(mol_act).readlines())
    all_inactive_mol = len(open(mol_inact).readlines())
    # ts_inactive_mol = len(open(ts_inact).readlines())

    model = os.path.basename(in_pma)
    with open(in_pma) as fpma:
        d = json.loads(fpma.readline().strip())
        labels = ''.join(i[0] for i in d['feature_coords'])

    with open(ts_act) as f:
        for x in f:
            ts_active_mol.append(x.strip().split('\t')[0])
        
    with open(in_act_screen) as fs:
        for column in fs:
            act_screen.append(column.strip().split('\t')[0])

    act_screen = set(act_screen).difference(ts_active_mol)

    with open(ts_inact) as f:
        for x in f:
            ts_inactive_mol.append(x.strip().split('\t')[0])

    with open(in_inact_screen) as fs:
        for column in fs:
            inact_screen.append(column.strip().split('\t')[0])

    inact_screen = set(inact_screen).difference(ts_inactive_mol)

    TP = len(act_screen)
    FP = len(inact_screen)
    FN = all_active_mol - len(ts_active_mol) - TP
    TN = all_inactive_mol - len(ts_inactive_mol) - FP

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
                labels]


def main(mol_act, mol_inact, ts_act, ts_inact, path_to_pma, path_to_screen, out_external):
    
    start_time = time.time()

    df_result = pd.DataFrame(columns=['model', 'TP', 'FP', 'precision', 'FPR', 'recall', 'F1', 'F2', 'F05', 'BA', 'EF', 'features'])


    for enum, in_pma in enumerate(sorted(os.listdir(path_to_pma))):
        ppath = (os.path.abspath(os.path.join(os.path.abspath(path_to_screen), 'screen_active_{}.txt'.format(in_pma.split('.')[0]))),
                 os.path.abspath(os.path.join(os.path.abspath(path_to_screen), 'screen_inactive_{}.txt'.format(in_pma.split('.')[0]))),
                 os.path.abspath(os.path.join(path_to_pma, in_pma)))

        result = get_external_stat(mol_act, mol_inact, ts_act, ts_inact, ppath[2], ppath[0], ppath[1])
        if result:
            df_result.loc[enum] = result
        else:
            continue

    df_result = df_result.sort_values(by=['recall', 'F05', 'F2'], ascending=False)

    df_result.to_csv(out_external, index=None, sep='\t')
    sys.stderr.write('{}: ({}s)\n\n'.format(os.path.basename(out_external), round(time.time()-start_time,3)))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-mol_act', '--input_active_mol', metavar='active.smi', required=True,
                        help='.smi file with active molecules.')
    parser.add_argument('-mol_inact', '--input_inactive_mol', metavar='inactive.smi', required=True,
                        help='.smi file with inactive molecules.')
    # parser.add_argument('-idb', '--in_active_database', metavar='active.db', required=True,
    #                     help='input DB SQLite file with activity moleculs')
    parser.add_argument('-ts_act', '--trainset_active_mol', metavar='active.smi', required=True,
                        help='txt file with active molecules from training set.')
    parser.add_argument('-ts_inact', '--trainset_inactive_mol', metavar='inactive.smi', required=True,
                        help='txt file with inactive molecules from training set.')
    parser.add_argument('-ppma', '--path_to_pma', metavar='models/pYY/', required=True,
                        help='path to pma files')
    parser.add_argument('-pscreen', '--path_to_screen', metavar='screen/pYY/', required=True,
                        help='path to screen')
    parser.add_argument('-o', '--out_external', metavar='external_subsetXX_pYY.txt', default=None,
                        help='output text file, which will contain: model, n_act, n_inact, precision, recall, '
                             'F1, BA, EF, features, n_act_train, n_inact_train. ')

    args = vars(parser.parse_args())
    for o, v in args.items():
        if o == "input_active_mol": mol_act = v
        if o == "input_inactive_mol": mol_inact = v
        # if o == "in_active_database": in_adb = v
        if o == "trainset_active_mol": ts_act = v
        if o == "trainset_inactive_mol": ts_inact = v
        if o == "path_to_pma": path_to_pma = v
        if o == "path_to_screen": path_to_screen = v
        if o == "out_external": out_external = v

    if out_external is None:
        out_external = os.path.join(os.path.split(os.path.dirname(mol_act))[0], 'result.txt')


    main(mol_act, mol_inact, ts_act, ts_inact, path_to_pma, path_to_screen, out_external)
    
