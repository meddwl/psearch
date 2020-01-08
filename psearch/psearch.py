#!/usr/bin/env python3
from __future__ import absolute_import
import os
import sys
import argparse
from multiprocessing import Pool

from scripts.screen_db import screen_db
from scripts.external_statistics import calc_stat
from scripts.gen_pharm_models import gen_pharm_models
from scripts.select_training_set_rdkit import trainingset_formation


def create_parser():
    parser = argparse.ArgumentParser(description='Pharmacophore model building', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-ma', '--active_mol', metavar='active.smi', required=True,
                        help='.smi file with active molecules.')
    parser.add_argument('-mi', '--inactive_mol', metavar='inactive.smi', required=True,
                        help='.smi file with inactive molecules.')
    parser.add_argument('-adb', '--input_active_db', metavar='input.db', required=True,
                        help='SQLite DB with active pharmacophores (feature coordinates).')
    parser.add_argument('-idb', '--input_inactive_db', metavar='input.db', required=True,
                        help='SQLite DB with inactive pharmacophores (feature coordinates).')
    parser.add_argument('-ts', '--mode_train_set', metavar='1 2', nargs='+', type=int, default=[1, 2],
                        help='1 - form a training set by Strategy 1 (a single training set from centroids of '
                             'individual clusters), '
                             '2 - form a training set by Strategy 2 (separate training set per each cluster), '
                             '1 2 - form a training sets by Strategy 1 and Strategy 2,')
    parser.add_argument('-pts', '--path_trainset', metavar='path/training/set', nargs='+', type=str, default=None,
                        help='If None, the path will be generated automatically. ')
    parser.add_argument('-pma', '--path_pma', metavar='path/pma/files', default=None,
                        help='If None, the path will be generated automatically. ')
    parser.add_argument('-ps', '--path_screen', metavar='path/pma/files', default=None,
                        help='If None, the path will be generated automatically. ')
    parser.add_argument('-u', '--upper', metavar='6', type=int, default=1000,
                        help='upper number of features used for generation of subpharmacophores.'
                             'if None well generated pharmacophores with as many features as possible')
    parser.add_argument('-tol', '--tolerance', metavar='VALUE', default=0,
                        help='tolerance used for calculation of a stereoconfiguration sign.')
    parser.add_argument('-thr', '--threshold_clust', default=0.4,
                        help='threshold for сlustering data by Butina algorithm')
    parser.add_argument('-f', '--min_features', metavar='INTEGER', default=None, type=int,
                        help='minimum number of features with distinct coordinates in models. '
                             'Default: all models will be screened.')
    parser.add_argument('--fdef', metavar='smarts.fdef', default=None,
                        help='fdef-file with pharmacophore feature definition.')
    parser.add_argument('-c', '--ncpu', metavar='cpu_number', default=1,
                        help='number of cpus to use for calculation.')
    return parser


def creating_pharmacophore(in_adb, in_indb, files_ats, files_ints, path_pma, tol, upper):
    path_pma, lower = gen_pharm_models(in_adb=in_adb, in_indb=in_indb,
                                       act_trainset=files_ats, inact_trainset=files_ints,
                                       out_pma=path_pma,
                                       tolerance=tol,
                                       lower=3, upper=upper)
    if lower == 0:
        sys.stderr.write('Error: {}.\nFor {},\n{}\n\n'.format(path_pma, files_ats, files_ints))
    if not path_pma:
        sys.stderr.write("Error: no folder with pma files.\nFor {},\n{}\n\n".format(files_ats, files_ints))


def pharmacophore_validation(mol_act, mol_inact, in_adb, in_indb, path_ts, path_pma, path_screen, min_features, ncpu):
    screen_act = os.path.join(path_screen, 'active')
    if not os.path.exists(screen_act):
        os.makedirs(screen_act)

    screen_inact = os.path.join(path_screen, 'inactive')
    if not os.path.exists(screen_inact):
        os.makedirs(screen_inact)

    screen_db(db_fname=in_adb, queries=path_pma, output=screen_act, input_sdf=None,
              match_first_conf=True, min_features=min_features, ncpu=ncpu)
    screen_db(db_fname=in_indb, queries=path_pma, output=screen_inact, input_sdf=None,
              match_first_conf=True, min_features=min_features, ncpu=ncpu)

    out_external_stat = os.path.join(os.path.split(os.path.dirname(os.path.abspath(in_adb)))[0], 'results',
                                     'external_statistics.txt')
    if not os.path.exists(os.path.dirname(out_external_stat)):
        os.makedirs(os.path.dirname(out_external_stat))

    calc_stat(mol_act=mol_act,
              mol_inact=mol_inact,
              path_ts=path_ts,
              path_to_pma=path_pma,
              in_act_screen=screen_act,
              in_inact_screen=screen_inact,
              out_external=out_external_stat)


def creating_pharmacophore_mp(items):
    new = creating_pharmacophore
    return creating_pharmacophore(*items)
    
    
def get_items(in_adb, in_indb, list_ts, path_pma, tol, upper):
    for file_ats, file_ints in list_ts:
        yield in_adb, in_indb, file_ats, file_ints, path_pma, tol, upper


def main(mol_act, mol_inact, in_adb, in_indb, mode_train_set, path_ts, path_pma, path_screen,
         tol, upper, min_features, fdef_fname, threshold_clust, ncpu):

    p = Pool(ncpu)

    if 1 in mode_train_set:
        list_ts_1 = trainingset_formation(in_fname_act=mol_act,
                                          in_fname_inact=mol_inact,
                                          output=path_ts,
                                          fdef_fname=fdef_fname,
                                          make_clust=False,
                                          fcfp4=False,
                                          clust_stat=None,
                                          threshold_clust=threshold_clust,
                                          clust_size=5,
                                          max_nact_trainset=5)
    else:
        list_ts_1 = []

    if 2 in mode_train_set:
        list_ts_2 = trainingset_formation(in_fname_act=mol_act,
                                          in_fname_inact=mol_inact,
                                          output=path_ts,
                                          fdef_fname=fdef_fname,
                                          make_clust=True,
                                          fcfp4=False,
                                          clust_stat=open(os.path.join(os.path.split(os.path.dirname(mol_act))[0],
                                               'cluster_stat_thr{}.txt'.format(threshold_clust)), 'w'),
                                          threshold_clust=threshold_clust,
                                          clust_size=5,
                                          max_nact_trainset=5)
    else:
        list_ts_2 = []

    list_ts = list_ts_1 + list_ts_2

    for _ in p.imap(creating_pharmacophore_mp, get_items(in_adb=in_adb, in_indb=in_indb, list_ts=list_ts,
                                                         path_pma=path_pma, tol=tol, upper=upper)):
        continue

    pharmacophore_validation(mol_act=mol_act,
                             mol_inact=mol_inact,
                             in_adb=in_adb,
                             in_indb=in_indb,
                             path_ts=path_ts,
                             path_pma=path_pma,
                             path_screen=path_screen,
                             min_features=min_features,
                             ncpu=ncpu)


def entry_point():
    parser = create_parser()
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
        if o == "upper": upper = int(v)
        if o == "threshold_clust": threshold_clust = float(v)
        if o == "min_features": min_features = int(v) if v is not None else None
        if o == "fdef": fdef_fname = v
        if o == "ncpu": ncpu = int(v)

    # creating paths for training sets, pma files and screening files
    if not path_ts:
        path_ts = os.path.join(os.path.split(os.path.dirname(os.path.abspath(in_adb)))[0], 'trainset')
        if not os.path.isdir(path_ts):
            os.makedirs(path_ts)

    if not path_pma:
        path_pma = os.path.join(os.path.split(os.path.dirname(in_adb))[0], 'models')
        if not os.path.isdir(path_pma):
            os.makedirs(path_pma)

    if not path_screen:
        path_screen = os.path.join(os.path.split(path_pma)[0], 'screen')
        if not os.path.isdir(path_screen):
            os.makedirs(path_screen)

    main(mol_act=active_mol,
         mol_inact=inactive_mol,
         in_adb=in_adb,
         in_indb=in_indb,
         path_ts=path_ts,
         path_pma=path_pma,
         path_screen=path_screen,
         mode_train_set=mode_train_set,
         tol=tol,
         upper=upper,
         threshold_clust=threshold_clust,
         min_features=min_features,
         fdef_fname=fdef_fname,
         ncpu=ncpu)


if __name__ == '__main__':
    entry_point()
