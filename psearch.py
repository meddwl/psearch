#!/usr/bin/env python3

import os
import time
import argparse

from scripts import gen_pharm_models
from scripts import select_training_set_rdkit
from scripts import screen_db
import external_statistics


def calc(mol_act, mol_inact, in_adb, in_indb, files_at, files_int, path_pma, path_screen, tol, lower):

    if len(os.listdir(path_pma)) == 0:
        path_pma, ids = gen_pharm_models.main(in_adb=in_adb,
                                              in_indb=in_indb,
                                              act_trainset=files_at,
                                              inact_trainset=files_int,
                                              out_pma=path_pma,
                                              tolerance=tol,
                                              lower=lower)
        if ids == 0:
            raise Exception('Error: {}'.format(path_pma))
        if not path_pma:
            raise Exception("Error: no folder with pma files")

    if len(os.listdir(path_pma)) != 0:
        start = time.time()
        path_screen = os.path.join(path_screen,
                                   '{}_ph{}'.format(os.path.basename(files_at).split('.')[0].split('_')[1], ids))
        if not os.path.exists(path_screen):
            os.mkdir(path_screen)

        screen_db.main(dbs_fname=[in_adb, in_indb], path_pma=path_pma, path_screen=path_screen)

        print('screen time : {}'.format(time.time() - start))

        out_external_stat = os.path.join(os.path.split(os.path.dirname(os.path.abspath(in_adb)))[0],
                                         'external_statistics{}.txt'.format(os.path.split(path_pma)[1]))

        external_statistics.main(mol_act=mol_act,
                                 mol_inact=mol_inact,
                                 ts_act=files_at,
                                 ts_inact=files_int,
                                 path_to_pma=path_pma,
                                 path_to_screen=path_screen,
                                 out_external=out_external_stat)


def main(mol_act, mol_inact, in_adb, in_indb, mode_train_set, path_ts, path_pma, path_screen,
         tol, lower, fdef_fname, threshold_clust, clust_size, max_nact_trainset):

    if type(path_ts) is str:
        if 1 in mode_train_set:
            list_ts_1 = select_training_set_rdkit.main(in_fname_act=mol_act,
                                                       in_fname_inact=mol_inact,
                                                       output=path_ts,
                                                       fdef_fname=fdef_fname,
                                                       make_clust=False,
                                                       fcfp4=fcfp4,
                                                       clust_stat=None,
                                                       threshold_clust=threshold_clust,
                                                       clust_size=clust_size,
                                                       max_nact_trainset=max_nact_trainset)
        else:
            list_ts_1 = []

        if 2 in mode_train_set:
            list_ts_2 = select_training_set_rdkit.main(in_fname_act=mol_act,
                                                       in_fname_inact=mol_inact,
                                                       output=path_ts,
                                                       fdef_fname=fdef_fname,
                                                       make_clust=True,
                                                       fcfp4=fcfp4,
                                                       clust_stat=open(os.path.join(os.path.split(os.path.dirname(mol_act))[0],
                                                           'cluster_stat_thr{}.txt'.format(threshold_clust)), 'w'),
                                                       threshold_clust=threshold_clust,
                                                       clust_size=clust_size,
                                                       max_nact_trainset=max_nact_trainset)
        else:
            list_ts_2 = []

        list_ts = list_ts_1 + list_ts_2

    else:
        list_ts = path_ts

    for files_at, files_int in list_ts:
        calc(mol_act, mol_inact,
             in_adb, in_indb,
             files_at, files_int,
             path_pma, path_screen,
             tol, lower)
        continue


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-ma', '--active_mol', metavar='active.smi', required=True,
                        help='.smi file with active molecules.')
    parser.add_argument('-mi', '--inactive_mol', metavar='inactive.smi', required=True,
                        help='.smi file with inactive molecules.')
    parser.add_argument('-adb', '--input_active_db', metavar='input.db', required=True,
                        help='SQLite DB with active pharmacophores (feature coordinates).')
    parser.add_argument('-idb', '--input_inactive_db', metavar='input.db', required=True,
                        help='SQLite DB with inactive pharmacophores (feature coordinates).')
    parser.add_argument('-ts', '--mode_train_set', nargs='+', default=[1, 2], type=int,
                        help='1 - form a training set by Stategy 1,'
                             '2 - form a training set by Stategy 2,'
                             '1 2 - form a training sets by Stategy 1 and Stategy 2,')
    parser.add_argument('-pts', '--path_trainset', metavar='path/training/set', nargs='+', type=str, default=None,
                        help='If None, the path will be generated automatically. '
                             'If the path to the folder is specified, and the folder is empty, '
                             'files with training set will be saved to the specified path. '
                             'or you can specify the path to the files with training sets, '
                             'where the first is a training set with active molecules, the second is with inactive ones.')
    parser.add_argument('-pma', '--path_pma', metavar='path/pma/files', default=None,
                        help='If None, the path will be generated automatically. '
                             'If the path to the folder is specified, and the folder is empty, '
                             'pma files (pharmacophore models) will be saved to the specified path. '
                             'If the path to the folder with the pma files is specified, '
                             'then their virtual screening will be conducted.')
    parser.add_argument('-ps', '--path_screen', metavar='path/pma/files', default=None,
                        help='If None, the path will be generated automatically. '
                             'If the path to the folder is specified, then screening files will be save in there')
    parser.add_argument('-tol', '--tolerance', metavar='VALUE', default=0,
                        help='tolerance used for calculation of a stereoconfiguration sign.')
    parser.add_argument('-l', '--lower', metavar='3', type=int, default=3,
                        help='lower number of features used for generation of subpharmacophores.')
    parser.add_argument('--fcfp4', action='store_true', default=False,
                        help='if set FCFP4 fingerprints will be used for compound selection, '
                             'otherwise pharmacophore fingerprints will be used based on feature '
                             'definitions provided by --rdkit_fdef argument.')
    parser.add_argument('-thr', '--threshold_clust', default=0.4,
                        help='threshold for сlustering data by Butina algorithm')
    parser.add_argument('-clz', '--clust_size', default=5,
                        help='minimum cluster size from extract centroinds for training set')
    parser.add_argument('-m', '--max_act_ts', default=5,
                        help='maximum number of active compounds for training set')
    parser.add_argument('--fdef', metavar='smarts.fdef',
                        default=os.path.join(os.getcwd(), 'pmapper/smarts_features.fdef'),
                        help='fdef-file with pharmacophore feature definition.')

    args = vars(parser.parse_args())
    for o, v in args.items():
        if o == "active_mol": active_mol = v
        if o == "inactive_mol": inactive_mol = v
        if o == "input_active_db": in_adb = v
        if o == "input_inactive_db": in_indb = v
        if o == "path_trainset": path_ts = v
        if o == "path_pma": path_pma = v
        if o == "path_screen": path_screen = v
        if o == "mode_train_set": mode_train_set = v
        if o == "tolerance": tol = int(v)
        if o == "lower": lower = int(v)
        if o == "input_active_mol": mol_act = v
        if o == "input_inactive_mol": mol_inact = v
        if o == "fcfp4": fcfp4 = v
        if o == "threshold_clust": threshold_clust = float(v)
        if o == "clust_size": clust_size = int(v)
        if o == "max_act_ts": max_nact_trainset = int(v)
        if o == "fdef": fdef_fname = v

    # creating paths for training sets, pma files and screening files
    if not path_ts:
        path_ts = os.path.join(os.path.split(os.path.dirname(os.path.abspath(in_adb)))[0], 'trainset')
        if not os.path.isdir(path_ts):
            os.makedirs(path_ts)
    elif len(path_ts) == 1 and os.path.isdir(path_ts[0]):
        path_ts = path_ts[0]
    elif os.path.isfile(path_ts[0]):
        path_ts = [[path_ts[i], path_ts[i+1]] for i in range(0, len(path_ts), 2)]

    if not path_pma:
        path_pma = os.path.join(os.path.split(os.path.dirname(in_adb))[0], 'models')
        if not os.path.isdir(path_pma):
            os.makedirs(path_pma)
    elif os.path.exists(path_pma) and os.path.isdir(path_pma):
        path_pma = path_pma

    if not path_screen:
        path_screen = os.path.join(os.path.split(path_pma)[0], 'screen')
        if not os.path.isdir(path_screen):
            os.makedirs(path_screen)
    elif os.path.exists(path_screen) and os.path.isdir(path_screen):
        path_screen = path_screen

    main(mol_act=active_mol,
         mol_inact=inactive_mol,
         in_adb=in_adb,
         in_indb=in_indb,
         path_ts=path_ts,
         path_pma=path_pma,
         path_screen=path_screen,
         mode_train_set=mode_train_set,
         tol=tol,
         lower=lower,
         threshold_clust=threshold_clust,
         clust_size=clust_size,
         max_nact_trainset=max_nact_trainset,
         fdef_fname=fdef_fname)