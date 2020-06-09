#!/usr/bin/env python3
# author          : Alina Kutlushina
# date            : 01.05.2018
# license         : BSD-3
#==============================================================================

import os
import argparse
import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.ML.Cluster import Butina
from rdkit.Chem import AllChem
from rdkit.Chem.Pharm2D import Generate
from rdkit.Chem.Pharm2D.SigFactory import SigFactory
from pmapper.customize import load_factory
from psearch.database import DB


def create_parser():
    parser = argparse.ArgumentParser(description='select compounds for training set',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--input_mols', metavar='input_molecules.smi', required=True,
                        help='path to input SMILES file with compounds.')
    parser.add_argument('-db', '--input_db', metavar='input_database.dat', required=True,
                        help='path to input database file.')
    parser.add_argument('-o', '--output', metavar='output/path', default=None,
                        help='output path. The folder where will be saved a training set.')
    parser.add_argument('-ts', '--mode_train_set', metavar='1 2', nargs='+', type=int, default=[1, 2],
                        help='Take numbers 1 or 2 or both to designate the strategy to create training sets. '
                             '1 - a single training set will be created from centroids of individual clusters, '
                             '2 - multiple training sets will be created, one per cluster. Default: 1 2.')
    parser.add_argument('--fcfp4', action='store_true', default=False,
                        help='if set FCFP4 fingerprints will be used for compound clustering, '
                             'otherwise pharmacophore fingerprints will be used.')
    parser.add_argument('-s', '--cluster_stat', default=None,
                        help='if designate path to file then save cluster statistics')
    parser.add_argument('-t', '--threshold_clust', type=float, default=0.4,
                        help='threshold for clustering data by Butina algorithm')
    parser.add_argument('-clz', '--clust_size', type=int, default=5,
                        help='minimum cluster size to extract centroids for the training set')
    parser.add_argument('-m', '--max_acts', type=int, default=5,
                        help='maximum number of active compounds for training set')
    return parser


def read_file(fname, db_fname, fcfp4):
    """
    :param fname: path to input SMILES file with molecules defined activity
    :param fcfp4:
    :return: pandas.DataFrame, columns = mol_name, smiles, activity, fp
    """
    df = pd.read_csv(fname, sep='\t', header=None)
    df.rename(columns={0: 'smiles', 1: 'mol_name', 2: 'activity'}, inplace=True)
    # drop columns if it is
    if not Chem.MolFromSmiles(df.at[0, 'smiles']):
        df.drop(index=0)
    db = DB(db_fname)
    mol_names = db.get_mol_names()
    df = df[df['mol_name'].isin(mol_names)]

    if fcfp4:
        df['fp'] = [(AllChem.GetMorganFingerprint(Chem.MolFromSmiles(smiles), 2, useFeatures=True)) for smiles in df['smiles']]
    else:
        featfactory = load_factory()
        sigfactory = SigFactory(featfactory, minPointCount=2, maxPointCount=3, trianglePruneBins=False)
        sigfactory.SetBins([(0, 2), (2, 5), (5, 8)])
        sigfactory.Init()
        df['fp'] = [(Generate.Gen2DFingerprint(Chem.MolFromSmiles(smiles), sigfactory)) for smiles in df['smiles']]
    return df


def gen_cluster_subset_butina(fps, cutoff):
    dists = []
    nfps = len(fps)
    for i in range(1, nfps):
        sims = DataStructs.BulkTanimotoSimilarity(fps[i], fps[:i])
        dists.extend([1 - x for x in sims])
    cs = Butina.ClusterData(dists, nfps, cutoff, isDistData=True)
    return cs  # returns tuple of tuples with sequential numbers of compounds in each cluster


def save_cluster_stat(cs_index, len_act, clust_stat):
    for i, cluster in enumerate(cs_index):
        i_act = 0
        for el in cluster:
            if el in range(len_act):
                i_act += 1
        # print('cluster №%i, cluster length %i, share of active %.2f' % (i, len(cluster), i_act/len(cluster)))
        # print(cluster, '\n')
        clust_stat.write(f'cluster {i}, cluster length {len(cluster)}, share of active {i_act/len(cluster)} \n')


def diff_binding_mode(cs, df_mols, len_act, inact_centroids, min_num):
    # mol_names contains actives and then inactives
    # therefore len_act is equal to the number of actives to select them from the list of mol names
    ts_full = []
    for i, c in enumerate(cs):
        if len(c) >= min_num:
            ts_mol_name_act = []
            ts_mol_name_inact = []
            for x in c:
                if x in range(len_act) and len(ts_mol_name_act) < min_num:
                    ts_mol_name_act.append((df_mols.at[x, 'mol_name'], df_mols.at[x, 'smiles']))
                elif x not in range(len_act) and len(ts_mol_name_inact) < min_num:
                    ts_mol_name_inact.append((df_mols.at[x, 'mol_name'], df_mols.at[x, 'smiles']))

            ts_mol_name_inact = set(ts_mol_name_inact) | set(inact_centroids)
            if len(ts_mol_name_act) == min_num:
                ts_full.append((i, tuple(sorted(ts_mol_name_act)), tuple(sorted(ts_mol_name_inact))))
    return tuple(ts_full)


def get_centroids(cs, df, num):
    return tuple(sorted((df.at[x[0], 'mol_name'], df.at[x[0], 'smiles']) for x in cs if len(x) >= num))


def trainingset_formation(input_mols, input_db, path_ts, mode_train_set, fcfp4,
                          clust_stat, threshold, clust_size, max_num_acts):

    if (1 not in mode_train_set) and (2 not in mode_train_set):
        return 'Wrong value of parameter mode_train_set. That should be 1 and/or 2.'

    df_mols = read_file(input_mols, input_db, fcfp4)
    if df_mols['activity'].dtypes == 'int64':
        df_mols = df_mols.sort_values(by='activity', ascending=True).reset_index(drop=True)
        inact_mark = 0
    else:
        df_mols = df_mols.sort_values(by='activity').reset_index(drop=True)
        inact_mark = 'inactive'
    len_acts = df_mols[df_mols['activity'] != inact_mark].shape[0]
    list_ts = []

    if 2 in mode_train_set:
        cs = gen_cluster_subset_butina(df_mols['fp'].tolist(), threshold)
        df_inact = df_mols[df_mols['activity'] == inact_mark].reset_index(drop=True)
        cs_inact = gen_cluster_subset_butina(df_inact['fp'].tolist(), threshold)
        inact_centroids = get_centroids(cs_inact, df_inact, clust_size)  # tuple of tuples with mol names and their SMILES

        if clust_stat:
            save_cluster_stat(cs, len_acts-1, clust_stat)

        ts_full = diff_binding_mode(cs, df_mols, len_acts-1, inact_centroids, max_num_acts)
        for i, act_ts, inact_ts in ts_full:
            out_act = os.path.join(path_ts, f'active_t{i}.smi')
            out_inact = os.path.join(path_ts, f'inactive_t{i}.smi')
            list_ts.append([out_act, out_inact])
            with open(out_act, 'wt') as f:
                f.write('\n'.join(f'{smiles}\t{mol_name}' for mol_name, smiles in act_ts))
            with open(out_inact, 'wt') as f:
                f.write('\n'.join(f'{smiles}\t{mol_name}' for mol_name, smiles in inact_ts))

    if 1 in mode_train_set:
        out_act = os.path.join(path_ts, 'active_centroid.smi')
        out_inact = os.path.join(path_ts, 'inactive_centroid.smi')

        # process actives
        df_act = df_mols[df_mols['activity'] != inact_mark].reset_index(drop=True)
        cs = gen_cluster_subset_butina(df_act['fp'].tolist(), threshold)
        centroids = get_centroids(cs, df_act, clust_size)
        if len(centroids) < clust_size:
            return list_ts
        with open(out_act, 'wt') as f:
            f.write('\n'.join(f'{smiles}\t{mol_name}' for mol_name, smiles in centroids))

        # process inactives
        df_inact = df_mols[df_mols['activity'] == inact_mark].reset_index(drop=True)
        cs = gen_cluster_subset_butina(df_inact['fp'].tolist(), threshold)
        centroids = get_centroids(cs, df_inact, clust_size)
        with open(out_inact, 'wt') as f:
            f.write('\n'.join(f'{smiles}\t{mol_name}' for mol_name, smiles in centroids))
        list_ts.append([out_act, out_inact])
    return list_ts


def entry_point():
    parser = create_parser()
    args = parser.parse_args()

    if not args.output:
        output = os.path.join(os.path.dirname(os.path.abspath(args.input_mols)), 'trainset')
    else:
        output = args.output
    if not os.path.exists(output):
        os.makedirs(output)

    trainingset_formation(input_mols=args.input_mols,
                          input_db=args.input_db,
                          path_ts=output,
                          mode_train_set=args.mode_train_set,
                          fcfp4=args.fcfp4,
                          clust_stat=open(args.cluster_stat, 'wt') if args.cluster_stat else None,
                          threshold=args.threshold_clust,
                          clust_size=args.clust_size,
                          max_num_acts=args.max_acts)


if __name__ == '__main__':
    entry_point()
