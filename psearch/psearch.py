#!/usr/bin/env python3
# author          : Alina Kutlushina
# date            : 10.01.2019
# license         : BSD-3
# ==============================================================================

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
    parser.add_argument('-p', '--project_dir', metavar='DIRNAME', type=str, default=None,
                        help='A path to a project dir. Directory where all intermediate and output files will be saved.')
    parser.add_argument('-i', '--molecules', metavar='molecules.smi', type=str, required=True,
                        help='The script takes as input a tab-separated SMILES file containing `SMILES`, `compound id`, '
                             '`activity` columns'
                             'The third column should contain a word 1 or 0. 1 is for actives, 0 is for inactives.')
    parser.add_argument('-d', '--database', metavar='FILENAME.dat', type=str, required=True,
                        help='Path to the database with precomputed conformers and pharmacophores for the same input file.')
    parser.add_argument('-ts', '--trainset', metavar='DIRNAME', type=str, default=None,
                        help='A path to the folder where will be saved a training set.'
                             'If omitted, the path will be generated automatically relative to project directory.')
    parser.add_argument('-q', '--query', metavar='DIRNAME', type=str, default=None,
                        help='A path to a folder where will be saved the created pharmacophore models.'
                             'If omitted, the path will be generated automatically relative to project directory.')
    parser.add_argument('-s', '--screening', metavar='DIRNAME', type=str, default=None,
                        help='A text (.txt) file which will store names of compounds which fit the model. In the case '
                             'multiple query models or directories were supplied as input'
                             'this should be the path to a directory where output files will be created to store '
                             'screening results. If multiple directories were specified as input the corresponding '
                             'directories will be created in the output dir. Names of created  directories will be '
                             'taken from the bottom level of input directories, e.g. path/to/model/ will be stored in '
                             'output_dir/model. Beware, existed output files/directories will be overwritten.')
    parser.add_argument('-r', '--external_statistics', metavar='FILENAME', default=None,
                        help='An output text file where will be saved validation statistics'
                             'If omitted, the path will be generated automatically relative to project directory.')
    parser.add_argument('-m', '--mode_train_set', nargs='+', type=int, default=[1, 2],
                        help='Take numbers 1 or 2 or both to designate the strategy to create training sets. '
                             '1 - a single training set will be created from centroids of individual clusters, '
                             '2 - multiple training sets will be created, one per cluster.')
    parser.add_argument('--fcfp4', action='store_true', default=False,
                        help='If set FCFP4 fingerprints will be used for compound clustering, '
                             'otherwise pharmacophore fingerprints will be used.')
    parser.add_argument('-t', '--threshold', metavar='NUMERIC', type=float, default=0.4,
                        help='threshold for сlustering data by Butina algorithm.')
    parser.add_argument('-tol', '--tolerance', metavar='NUMERIC', type=float, default=0,
                        help='tolerance used for calculation of a stereoconfiguration sign.')
    parser.add_argument('-b', '--bin_step', metavar='NUMERIC', type=float, default=1,
                        help='binning step.')
    parser.add_argument('-l', '--lower', metavar='INTEGER', type=int, default=3,
                        help='starting from this number of features, pharmacophore models will be created')
    parser.add_argument('-f', '--save_model_complexity', metavar='INTEGER', type=int, default=None,
                        help='All pharmacophore models will be saved starting from this number of features.'
                             'If omitted will be saved only the most complex pharmacophore models')
    parser.add_argument('-u', '--upper', metavar='INTEGER', type=int, default=None,
                        help='limit the upper number of features in generated pharmacophores. '
                             'If omitted pharmacophores of maximum complexity will be generated.')
    parser.add_argument('-c', '--ncpu', metavar='INTEGER', type=int, default=1,
                        help='number of cpus to use for calculation.')
    return parser


def creating_pharmacophore_mp(items):
    return creating_pharmacophore(*items)


def get_items(in_db, list_ts, path_pma, upper, lower, save_model_complexity, bin_step, tolerance):
    for train_set in list_ts:
        yield in_db, train_set, path_pma, upper, lower, save_model_complexity, bin_step, tolerance


def creating_pharmacophore(in_db, train_set, path_pma, upper, lower, save_model_complexity, bin_step, tolerance):
    gen_pharm_models(in_db=in_db,
                     trainset=train_set,
                     out_pma=path_pma,
                     upper=upper,
                     bin_step=bin_step,
                     tolerance=tolerance,
                     current_nfeatures=lower,
                     nfeatures=save_model_complexity,
                     save_statistics=False)


def pharmacophore_validation(mols, in_db, path_ts, path_pma, path_screen, pp_external_stat):
    screen_db(db_fname=in_db, queries=[os.path.join(path_pma, mm) for mm in os.listdir(path_pma)],
              output=path_screen, output_sdf=None,
              match_first_conf=True, min_features=None, ncpu=1)

    calc_stat(path_mols=mols,
              path_ts=path_ts,
              pp_models=path_pma,
              path_screen=path_screen,
              out_external=pp_external_stat)


def main(in_mols, in_db, path_ts, path_pma, path_screen, path_external_stat, path_clus_stat,
         mode_train_set, fcfp4, threshold, tolerance, lower, save_model_complexity, upper, bin_step, ncpu):

    p = Pool(ncpu)
    list_ts = trainingset_formation(input_mols=in_mols,
                                    path_ts=path_ts,
                                    mode_train_set=mode_train_set,
                                    fcfp4=fcfp4,
                                    clust_stat=open(path_clus_stat, 'w'),
                                    threshold=threshold)

    if type(list_ts) == str:
        sys.exit(list_ts)

    for _ in p.imap(creating_pharmacophore_mp, get_items(in_db=in_db, list_ts=list_ts, path_pma=path_pma,
                                                         upper=upper, lower=lower,
                                                         save_model_complexity=save_model_complexity,
                                                         bin_step=bin_step, tolerance=tolerance)):
        continue

    pharmacophore_validation(mols=in_mols,
                             in_db=in_db,
                             path_ts=path_ts,
                             path_pma=path_pma,
                             path_screen=path_screen,
                             pp_external_stat=path_external_stat)


def entry_point():
    parser = create_parser()
    args = parser.parse_args()
    project_dir = args.project_dir if args.project_dir else os.path.dirname(os.path.abspath(args.molecules))
    os.makedirs(project_dir, exist_ok=True)
    main(in_mols=os.path.abspath(args.molecules),
         in_db=os.path.abspath(args.database),
         path_ts=args.trainset if args.trainset else os.path.join(project_dir, 'trainset'),
         path_pma=args.query if args.query else os.path.join(project_dir, 'models'),
         path_screen=args.screening if args.screening else os.path.join(project_dir, 'screen'),
         path_external_stat=args.external_statistics if args.external_statistics else os.path.join(project_dir, 'results', 'external_statistics.txt'),
         path_clus_stat=os.path.join(project_dir, 'cluster_stat_t{}.txt'.format(args.threshold)),
         mode_train_set=args.mode_train_set,
         fcfp4=args.fcfp4,
         threshold=float(args.threshold),
         tolerance=float(args.tolerance),
         bin_step=int(args.bin_step),
         lower=int(args.lower),
         save_model_complexity=int(args.save_model_complexity) if args.save_model_complexity else None,
         upper=int(args.upper) if args.upper is not None else None,
         ncpu=int(args.ncpu))


if __name__ == '__main__':
    entry_point()
