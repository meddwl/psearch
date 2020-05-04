#!/usr/bin/env python3
import os
import sys
import argparse
from multiprocessing import Pool

from psearch.screen_db import screen_db
from psearch.scripts.external_statistics import calc_stat
from psearch.scripts.gen_pharm_models import gen_pharm_models
from psearch.scripts.select_training_set_rdkit import trainingset_formation


def create_parser():
    parser = argparse.ArgumentParser(description='Ligand-based pharmacophore model building',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-p', '--project_dir', type=str, default=None,
                        help='path to a project dir. Path where will be saved ')
    parser.add_argument('-m', '--input_molecules', metavar='molecules.smi', type=str, required=True,
                        help='path to tab-separated file with SMILES, molecule name and active/inactive '
                             'in the third column.')
    parser.add_argument('-db', '--input_db', metavar='FILENAME', type=str, required=True,
                        help='path to shelve database with precomputed conformers and pharmacophores')
    parser.add_argument('-ts', '--mode_train_set', metavar='1 2', nargs='+', type=int, default=[1, 2],
                        help='Take numbers 1 or 2 or both to designate the strategy to create training sets. '
                             '1 - a single training set will be created from centroids of individual clusters, '
                             '2 - multiple training sets will be created, one per cluster. Default: 1 2.')
    parser.add_argument('-b', '--bin_step', metavar='NUMERIC', type=float, default=1,
                        help='binning step. Default: 1.')
    parser.add_argument('-u', '--upper', metavar='INTEGER', type=int, default=None,
                        help='upper number of features in generation of pharmacophores.'
                             'if omitted pharmacophores of maximum complexity will be generated..')
    parser.add_argument('-tol', '--tolerance', metavar='NUMERIC', type=float, default=0,
                        help='tolerance used for calculation of a stereoconfiguration sign.')
    parser.add_argument('-thr', '--threshold', metavar='NUMERIC', type=float, default=0.4,
                        help='threshold for сlustering data by Butina algorithm.')
    parser.add_argument('-pts', '--path_trainset', metavar='path/training/set', type=str, default=None,
                        help='If omitted, the path will be generated automatically.')
    parser.add_argument('-pm', '--path_models', metavar='path/to/models/', type=str, default=None,
                        help='If omitted, the path will be generated automatically.')
    parser.add_argument('-ps', '--path_screen', metavar='path/to/screen/output', type=str, default=None,
                        help='If omitted, the path will be generated automatically.')
    parser.add_argument('-pr', '--path_external_statistics', metavar='path/external/statistics', default=None,
                        help='If omitted, the path will be generated automatically.')
    parser.add_argument('-c', '--ncpu', metavar='cpu_number', default=1,
                        help='number of cpus to use for calculation.')
    return parser


def creating_pharmacophore_mp(items):
    return creating_pharmacophore(*items)


def get_items(project_dir, in_db, list_ts, path_pma, tol, upper, bin_step):
    for file_ats, file_ints in list_ts:
        yield project_dir, in_db, file_ats, file_ints, path_pma, tol, upper, bin_step


def creating_pharmacophore(project_dir, in_db, files_ats, files_ints, path_pma, tol, upper, bin_step):
    gen_pharm_models(project_dir=project_dir,
                     in_db=in_db,
                     act_trainset=files_ats,
                     inact_trainset=files_ints,
                     out_pma=path_pma,
                     tolerance=tol,
                     lower=3, upper=upper,
                     bin_step=bin_step,
                     save_models=0,
                     save_files=True)


def pharmacophore_validation(mols, in_db, path_ts, path_pma, path_screen, pp_external_stat):

    screen_db(db_fname=in_db, queries=[path_pma], output=path_screen, output_sdf=None,
              match_first_conf=True, min_features=None, ncpu=1)

    calc_stat(path_mols=mols,
              path_ts=path_ts,
              path_pma=path_pma,
              path_screen=path_screen,
              out_external=pp_external_stat)


def main(project_dir, in_mols, in_db, path_ts, path_pma, path_screen, pp_external_stat, path_clus_stat,
         mode_train_set, tol, upper, threshold, bin_step, ncpu):

    p = Pool(ncpu)

    if not os.path.exists(path_ts):
        os.makedirs(path_ts)
    list_ts = trainingset_formation(input_mols=in_mols,
                                    path_ts=path_ts,
                                    # fdef_fname=fdef_fname,
                                    mode_train_set=mode_train_set,
                                    fcfp4=False,
                                    clust_stat=open(path_clus_stat, 'w'),
                                    threshold=threshold,
                                    clust_size=5,
                                    max_num_acts=5)

    if type(list_ts) == str:
        sys.exit(list_ts)

    for _ in p.imap(creating_pharmacophore_mp, get_items(project_dir=project_dir, in_db=in_db, list_ts=list_ts,
                                                         path_pma=path_pma, tol=tol,
                                                         upper=upper, bin_step=bin_step)):
        continue

    pharmacophore_validation(mols=in_mols,
                             in_db=in_db,
                             path_ts=path_ts,
                             path_pma=path_pma,
                             path_screen=path_screen,
                             pp_external_stat=pp_external_stat)


def entry_point():
    parser = create_parser()
    args = vars(parser.parse_args())
    for o, v in args.items():
        if o == "project_dir": project_dir = v
        if o == "input_molecules": in_mols = os.path.abspath(v)
        if o == "input_db": in_db = os.path.abspath(v)
        if o == "path_trainset": path_ts = v
        if o == "path_models": path_pma = v
        if o == "path_screen": path_screen = v
        if o == "path_external_statistics": pp_external_stat = v
        if o == "mode_train_set": mode_train_set = v
        if o == "tolerance": tol = int(v)
        if o == "upper": upper = int(v) if v is not None else None
        if o == "threshold": threshold = float(v)
        # if o == "fdef": fdef_fname = v
        if o == "bin_step": bin_step = int(v)
        if o == "ncpu": ncpu = int(v)

    # creating paths for training sets, pma files and screening files
    if not project_dir:
        project_dir = os.path.dirname(in_mols)

    if not path_ts:
        path_ts = os.path.join(project_dir, 'trainset')

    if not path_pma:
        path_pma = os.path.join(project_dir, 'models')

    if not path_screen:
        path_screen = os.path.join(project_dir, 'screen')

    if not pp_external_stat:
        pp_external_stat = os.path.join(project_dir, 'results', 'external_statistics.txt')

    main(project_dir=project_dir,
         in_mols=in_mols,
         in_db=in_db,
         path_ts=path_ts,
         path_pma=path_pma,
         path_screen=path_screen,
         pp_external_stat=pp_external_stat,
         path_clus_stat=os.path.join(project_dir, 'cluster_stat_t{}.txt'.format(threshold)),
         mode_train_set=mode_train_set,
         tol=tol,
         upper=upper,
         threshold=threshold,
         bin_step=bin_step,
         ncpu=ncpu)


if __name__ == '__main__':
    entry_point()
