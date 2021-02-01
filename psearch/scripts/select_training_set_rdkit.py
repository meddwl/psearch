#!/usr/bin/env python3
# author          : Alina Kutlushina
# date            : 01.05.2018
# license         : BSD-3
#==============================================================================

import os
import argparse
import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.ML.Cluster import Butina
from rdkit.Chem import AllChem
from rdkit.Chem.Pharm2D import Generate
from rdkit.Chem.Pharm2D.SigFactory import SigFactory
from pmapper.customize import load_factory


def create_parser():
    parser = argparse.ArgumentParser(description='select compounds for training set',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--input_mols', metavar='input_molecules.smi', required=True,
                        help='The script takes as input a tab-separated SMILES file containing `SMILES`, '
                             '`compound id`, `activity` columns. '
                             'The third column should contain a word 1 or 0. 1 is for actives, 0 is for inactive ones.')
    parser.add_argument('-o', '--output', metavar='output/path', default=None,
                        help='An output path to the folder where will be saved a training set.'
                             'If omitted, the path will be generated automatically relative to project directory.')
    parser.add_argument('-ts', '--mode_train_set', metavar='1 2', nargs='+', type=int, default=[1, 2],
                        help='Take numbers 1 or 2 or both to designate the strategy to create training sets. '
                             '1 - a single training set will be created from centroids of individual clusters, '
                             '2 - multiple training sets will be created, one per cluster.')
    parser.add_argument('--fcfp4', action='store_true', default=False,
                        help='If set FCFP4 fingerprints will be used for compound clustering, '
                             'otherwise pharmacophore fingerprints will be used.')
    parser.add_argument('-thr', '--threshold_clust', type=float, default=0.4,
                        help='threshold for clustering data by Butina algorithm')
    parser.add_argument('-s', '--save_statistics', default=None,
                        help='If writen path to file then cluster statistics will be saved into this file')
    return parser


def read_file(fname, fcfp4):
    df = pd.read_csv(fname, sep='\t', header=None)
    df.rename(columns={0: 'smiles', 1: 'mol_name', 2: 'activity'}, inplace=True)
    if not Chem.MolFromSmiles(df.at[0, 'smiles']):
        df.drop(index=0, inplace=True)
    df['activity'] = df['activity'].astype(np.int64)
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


def save_cluster_stat(cs_index, df, index_acts, clust_stat):
    clust_stat.write("""#The file contains statistical information about the cluster 
#and the index of the molecule according to its location in the input file \n""")
    for i, cluster in enumerate(cs_index):
        i_act = 0
        for el in cluster:
            if el in index_acts:
                i_act += 1
        clust_stat.write(f'\ncluster {i}, cluster length {len(cluster)}, share of active {round(i_act/len(cluster), 3)}\n')
        clust_stat.write(','.join(map(str, [df.at[x, 'index'] for x in cluster])) + '\n')


def diff_binding_mode(cs, df_mols, index_acts, inact_centroids, min_num):
    ts_full = []
    for i, c in enumerate(cs):
        if len(set(c).intersection(index_acts)) >= min_num:
            ts_mol_name_act = tuple((df_mols.at[x, 'mol_name'], df_mols.at[x, 'smiles']) for x in list(set(c).intersection(index_acts))[:5])
            ts_mol_name_inact = [(df_mols.at[x, 'mol_name'], df_mols.at[x, 'smiles']) for x in list(set(c).difference(index_acts))[:5]]
            ts_mol_name_inact = tuple(set(ts_mol_name_inact).union(set(inact_centroids)))
            ts_full.append((i, ts_mol_name_act, ts_mol_name_inact))
    return tuple(ts_full)


def get_centroids(cs, df, num):
    return tuple(sorted((df.at[x[0], 'mol_name'], df.at[x[0], 'smiles']) for x in cs if len(x) >= num))


def trainingset_formation(input_mols, path_ts, mode_train_set, fcfp4, clust_stat, threshold):
    os.makedirs(path_ts, exist_ok=True)
    clust_size, max_num_acts = 5, 5
    if (1 not in mode_train_set) and (2 not in mode_train_set):
        return 'Wrong value of parameter mode_train_set. That should be 1 and/or 2.'

    df_mols = read_file(input_mols, fcfp4)
    df_inact = df_mols[df_mols['activity'] == 0].reset_index()
    df_act = df_mols[df_mols['activity'] == 1].reset_index()
    df_mols = df_act.append(df_inact).reset_index(drop=True)

    list_ts = []
    if 2 in mode_train_set:
        cs = gen_cluster_subset_butina(df_mols['fp'].tolist(), threshold)
        cs_inact = gen_cluster_subset_butina(df_inact['fp'].tolist(), threshold)
        centroids_inact = get_centroids(cs_inact, df_inact, clust_size)  # tuple of tuples with mol names and their SMILES

        if clust_stat:
            save_cluster_stat(cs, df_mols, range(df_act.shape[0]), clust_stat)

        ts_full = diff_binding_mode(cs, df_mols, range(df_act.shape[0]), centroids_inact, max_num_acts)
        for i, act_ts, inact_ts in ts_full:
            output = os.path.join(path_ts, f't{i}.smi')
            list_ts.append(output)
            with open(output, 'wt') as f:
                f.write('\n'.join(f'{smiles}\t{mol_name}\t{1}' for mol_name, smiles in act_ts) + '\n')
                f.write('\n'.join(f'{smiles}\t{mol_name}\t{0}' for mol_name, smiles in inact_ts))

    if 1 in mode_train_set:
        output = os.path.join(path_ts, 'centroids.smi')
        # process actives
        cs_act = gen_cluster_subset_butina(df_act['fp'].tolist(), threshold)
        centroids_act = get_centroids(cs_act, df_act, clust_size)
        if len(centroids_act) < max_num_acts:
            return list_ts
        # process inactives
        if 2 not in mode_train_set:
            cs_inact = gen_cluster_subset_butina(df_inact['fp'].tolist(), threshold)
            centroids_inact = get_centroids(cs_inact, df_inact, clust_size)
        with open(output, 'wt') as f:
            f.write('\n'.join(f'{smiles}\t{mol_name}\t{1}' for mol_name, smiles in centroids_act) + '\n')
            f.write('\n'.join(f'{smiles}\t{mol_name}\t{0}' for mol_name, smiles in centroids_inact))
        list_ts.append(output)
    return list_ts


def entry_point():
    parser = create_parser()
    args = parser.parse_args()
    trainingset_formation(input_mols=os.path.abspath(args.input_mols),
                          path_ts=args.output if args.output else os.path.join(os.path.dirname(os.path.abspath(args.input_mols)), 'trainset'),
                          mode_train_set=args.mode_train_set,
                          fcfp4=args.fcfp4,
                          clust_stat=open(args.save_statistics, 'wt') if args.save_statistics else None,
                          threshold=args.threshold_clust)


if __name__ == '__main__':
    entry_point()
