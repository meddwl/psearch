#!/usr/bin/env python3
# author          : Alina Kutlushina
# date            : 01.05.2018
# license         : BSD-3
#==============================================================================

import os
import argparse
import time
import scripts.create_subpharm as create_subpharm
from scripts import select_training_set_rdkit
import scripts.gen_subph as gen_subph
import scripts.screen_db as screen_db
import external_statistics
from multiprocessing import Pool


def calc(mol_act, mol_inact,
         in_adb, in_indb,
         files_ats, files_ints,
         tol, freq, lower):

    ph_act_subset_pX = create_subpharm.main(in_adb, files_ats, tol, lower, freq)
    ph_inact_subset_pX = create_subpharm.main(in_indb, files_ints, tol, lower, freq)

    path_pma, num = gen_subph.main(sub_act=ph_act_subset_pX,
                                   sub_inact=ph_inact_subset_pX,
                                   in_adb=in_adb,
                                   in_indb=in_indb,
                                   act_trainset=files_ats,
                                   lower=lower)

    if path_pma != None:
        
        path_screen = os.path.join(os.path.dirname(os.path.abspath(in_adb)), 'screen',
                                   os.path.basename(files_ats).split('.')[0].split('_')[1], 'pharm%i' % num)
        os.makedirs(path_screen, exist_ok=True)

        time_start = time.time()
        
        for query_fname in os.listdir(path_pma):
            out_f_screen_act = os.path.join(path_screen, 'screen_{}_{}.txt'.format(os.path.basename(in_adb).split('.')[0], query_fname.split('.')[0]))
            screen_db.main(db_fname=in_adb, query_fname=os.path.join(path_pma, query_fname), out_fname=out_f_screen_act, verbose=True, num=num)
            
            out_f_screen_inact = os.path.join(path_screen, 'screen_{}_{}.txt'.format(os.path.basename(in_indb).split('.')[0], query_fname.split('.')[0]))
            screen_db.main(db_fname=in_indb, query_fname=os.path.join(path_pma, query_fname), out_fname=out_f_screen_inact, verbose=True, num=num)
        print('screening of databases: ({}s)'.format(round(time.time()-time_start), 3))
        
        out_external_stat = os.path.join(os.path.dirname(os.path.abspath(in_adb)),
                                         'external_statistics{}_{}.txt'.format(
                                            os.path.basename(files_ats).split('.')[0].split('active')[1],
                                             os.path.split(path_pma)[1]))

        external_statistics.main(mol_act=mol_act,
                                 mol_inact=mol_inact,
                                 ts_act=files_ats,
                                 ts_inact=files_ints,
                                 path_to_pma=path_pma,
                                 path_to_screen=path_screen,
                                 out_external=out_external_stat)


def calc_mp(items):
    return calc(*items)
    
    
def get_items(mol_act, mol_inact,
             in_adb, in_indb,
             files_ats, files_ints,
             tol, freq, lower):
    for file_at, file_int in zip(files_ats, files_ints):
        yield mol_act, mol_inact, in_adb, in_indb, file_at, file_int, tol, freq, lower


def main(mol_act, mol_inact,
         in_adb, in_indb, mode_train_set,
         tol, freq,
         lower, fdef_fname, ncpu, treshold_clust, clust_size, max_nact_trainset):
    
    p = Pool(ncpu)

    if 1 in mode_train_set:
        list_ts_1 = select_training_set_rdkit.main(in_fname_act=mol_act,
                                                   in_fname_inact=mol_inact,
                                                   fdef_fname=fdef_fname,
                                                   make_clust=False,
                                                   fcfp4=fcfp4,
                                                   clust_stat=open(os.path.join(os.path.dirname(os.path.abspath(mol_act)), 'cluster_stat.txt'), 'wt'),
                                                   treshold_clust=treshold_clust,
                                                   clust_size=clust_size,
                                                   max_nact_trainset=max_nact_trainset)
    else:
        list_ts_1 = []

    if 2 in mode_train_set:
        list_ts_2 = select_training_set_rdkit.main(in_fname_act=mol_act,
                                                   in_fname_inact=mol_inact,
                                                   fdef_fname=fdef_fname,
                                                   make_clust=True,
                                                   fcfp4=fcfp4,
                                                   clust_stat=open(os.path.join(os.path.dirname(os.path.abspath(mol_act)),
                                                                    'cluster_stat_thr{}.txt'.format(treshold_clust)), 'wt'),
                                                   treshold_clust=treshold_clust,
                                                   clust_size=clust_size,
                                                   max_nact_trainset=max_nact_trainset)
    else:
        list_ts_2 = []

    list_ts = list_ts_1 + list_ts_2

    files_ats = [fts[0] for fts in list_ts]
    files_ints = [fts[1] for fts in list_ts]

    for _ in p.imap(calc_mp, get_items(mol_act, mol_inact,
                                       in_adb, in_indb,
                                       files_ats, files_ints,
                                       tol, freq,
                                       lower)):
        continue

    
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-p', '--project_dir', default=None,
                        help='path to a project dir.')
    parser.add_argument('-mol_act', '--active_mol', metavar='active.smi', default=None,
                        help='.smi file with active molecules. Can be omitted if project_dir is specified.')
    parser.add_argument('-mol_inact', '--inactive_mol', metavar='inactive.smi', default=None,
                        help='.smi file with inactive molecules. Can be omitted if project_dir is specified.')
    parser.add_argument('-iadb', '--input_active_db', metavar='input.db', default=None,
                        help='SQLite DB with active pharmacophores (feature coordinates). '
                             'Can be omitted if project_dir is specified.')
    parser.add_argument('-iindb', '--input_inactive_db', metavar='input.db', default=None,
                        help='SQLite DB with inactive pharmacophores (feature coordinates). '
                             'Can be omitted if project_dir is specified.')
    parser.add_argument('-ts', '--mode_train_set', nargs='+', default=[1, 2], type=int,
                        help='1 - form a training set by Stategy 1,'
                             '2 - form a training set by Stategy 2,'
                             '1 2 - form a training sets by Stategy 1 and Stategy 2,')
    parser.add_argument('-tol', '--tolerance', metavar='VALUE', default=0,
                        help='tolerance used for calculation of a stereoconfiguration sign.')
    parser.add_argument('-l', '--lower', metavar='4', type=int, default=4,
                        help='lower number of features used for generation of subpharmacophores.')
    parser.add_argument('-f', '--freq', metavar='0', default=None,
                        help='minimal frequency of pharmacophore occurrence in input compounds to keep them '
                             'in output file. Default: None (pharmacophores occurred in at least a half of '
                             'input compounds will be stored).')
    parser.add_argument('-c', '--ncpu', metavar='cpu_number', default=1,
                        help='number of cpus to use for calculation.')
    parser.add_argument('--fcfp4', action='store_true', default=False,
                        help='if set FCFP4 fingerprints will be used for compound selection, '
                             'otherwise pharmacophore fingerprints will be used based on feature '
                             'definitions provided by --rdkit_fdef argument.')
    parser.add_argument('-t', '--treshold_clust', default=0.4,
                        help='treshold for сlustering data by Butina algorithm')
    parser.add_argument('-cz', '--clust_size', default=5,
                        help='minimum cluster size from extract centroinds for training set')
    parser.add_argument('-ma', '--max_act_ts', default=5,
                        help='maximum number of active compounds for training set')
    parser.add_argument('--fdef', metavar='smarts.fdef', required=os.path.join(os.path.dirname(os.path.abspath(__file__)), 'pmapper', 'smarts_features.fdef'),
                        help='fdef-file with pharmacophore feature definition.')

    args = vars(parser.parse_args())
    for o, v in args.items():
        if o == "project_dir": project_dir = v
        if o == "active_mol": active_mol = v
        if o == "inactive_mol": inactive_mol = v
        if o == "input_active_db": in_adb = v
        if o == "input_inactive_db": in_indb = v
        if o == "mode_train_set": mode_train_set = v
        if o == "tolerance": tol = int(v)
        if o == "lower": lower = int(v)
        if o == "freq": freq = int(v) if v else None
        if o == "input_active_mol": mol_act = v
        if o == "input_inactive_mol": mol_inact = v
        if o == "fcfp4": fcfp4 = v
        if o == "ncpu": ncpu = int(v)
        if o == "treshold_clust": treshold_clust = float(v)
        if o == "clust_size": clust_size = int(v)
        if o == "max_act_ts": max_nact_trainset = int(v)
        if o == "fdef": fdef_fname = v

    if project_dir is not None:
        active_mol = os.path.join(project_dir, 'compounds', 'active.smi')
        inactive_mol = os.path.join(project_dir, 'compounds', 'inactive.smi')
        in_adb = os.path.join(project_dir, 'active.db')
        in_indb = os.path.join(project_dir, 'inactive.db')
    else:
        if active_mol is None or inactive_mol is None or in_adb is None or in_indb is None:
            print('If project dir was not set, active_mol, inactive_mol, input_active_db, input_inactive_db arguments '
                  'must be specified.')
            exit()

    main(mol_act=active_mol,
         mol_inact=inactive_mol,
         in_adb=in_adb,
         in_indb=in_indb,
         mode_train_set=mode_train_set,
         tol=tol,
         freq=freq,
         lower=lower,
         ncpu=ncpu,
         treshold_clust=treshold_clust,
         clust_size=clust_size,
         max_nact_trainset=max_nact_trainset,
         fdef_fname=fdef_fname)

