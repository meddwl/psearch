#!/usr/bin/env python3

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
    parser.add_argument('-i', '--input_mols', metavar='FILENAME.smi', required=True,
                        help='The script takes as input a tab-separated SMILES file containing `SMILES`, '
                             '`compound id`, `activity` columns. '
                             'The third column should contain a word 1 or 0. 1 is for actives, 0 is for inactive ones.')
    parser.add_argument('-o', '--output', metavar='DIRNAME', default=None,
                        help='An output path to the folder where will be saved a training set.'
                             'If omitted, the path will be generated automatically relative to project directory.')
    parser.add_argument('-ts', '--mode_train_set', metavar='1 2', nargs='+', type=int, default=[1, 2],
                        help='Take numbers 1 or 2 or both to designate the strategy to create training sets. '
                             '1 - a single training set will be created from centroids of individual clusters, '
                             '2 - multiple training sets will be created, one per cluster.')
    parser.add_argument('--fcfp4', action='store_true', default=False,
                        help='If set FCFP4 fingerprints will be used for compound clustering, '
                             'otherwise pharmacophore fingerprints will be used.')
    parser.add_argument('-t', '--threshold_clust', metavar='NUMERIC', type=float, default=0.4,
                        help='threshold for clustering data by Butina algorithm')
    parser.add_argument('-s', '--save_statistics', metavar='FILENAME', default=None,
                        help='If writen path to file then cluster statistics will be saved into this file')
    return parser


def read_file(fname, fcfp4):

    df = pd.read_csv(fname, sep='\t', header=None)
    df.rename(columns={0: 'smiles', 1: 'mol_name', 2: 'activity'}, inplace=True)
    if not Chem.MolFromSmiles(df.at[0, 'smiles']):
        df.drop(index=0, inplace=True)

    df['activity'] = df['activity'].astype(np.int64)
    df = df.sort_values(by=['activity'], ascending=False).reset_index(drop=True)
    df.to_csv(os.path.splitext(fname)[0] + '_sorted.smi', sep='\t', index=None)

    if fcfp4:
        fp = [(AllChem.GetMorganFingerprint(Chem.MolFromSmiles(smiles), 2, useFeatures=True)) for smiles in df['smiles']]
    else:
        featfactory = load_factory()
        sigfactory = SigFactory(featfactory, minPointCount=2, maxPointCount=3, trianglePruneBins=False)
        sigfactory.SetBins([(0, 2), (2, 5), (5, 8)])
        sigfactory.Init()
        fp = [(Generate.Gen2DFingerprint(Chem.MolFromSmiles(smiles), sigfactory)) for smiles in df['smiles']]
    return df, fp


def gen_cluster_subset_butina(fps, cutoff):
    dists = []
    for i in range(len(fps)-1):
        sims = DataStructs.BulkTanimotoSimilarity(fps[i], fps[i+1:])
        dists.extend([1 - x for x in sims])
    cs = Butina.ClusterData(dists, len(fps), cutoff, isDistData=True)
    return cs  # returns tuple of tuples with sequential numbers of compounds in each cluster


def save_cluster_stat(cs_index, df, index_acts, clust_stat):
    clust_stat.write("""#The saved file contains statistical information about the cluster 
#and the index of the molecule according to its location in the input file \n""")
    for i, cluster in enumerate(cs_index):
        i_act = len(set(cluster).intersection(index_acts))
        clust_stat.write(f'\ncluster {i}, cluster length {len(cluster)}, share of active {round(i_act/len(cluster), 3)}\n')
        clust_stat.write(','.join(map(str, cluster)) + '\n')


def save_train_set(output, act_ts, inact_ts):
    with open(output, 'wt') as f:
        f.write('\n'.join(f'{smiles}\t{mol_name}\t{activity}' for smiles, mol_name, activity in act_ts) + '\n')
        f.write('\n'.join(f'{smiles}\t{mol_name}\t{activity}' for smiles, mol_name, activity in inact_ts) + '\n')


def diff_binding_mode(cs, df, index_acts, inact_centroids, min_num):
    for i, c in enumerate(cs):
        if len(set(c).intersection(index_acts)) >= min_num:
            dfc = df.iloc[list(c)]
            ts_mol_name_act = dfc[dfc['activity'] == 1][:5].values
            ts_mol_name_inact = np.append(dfc[dfc['activity'] == 0][:5].values, inact_centroids, axis=0)
            ts_mol_name_inact = set(tuple(element) for element in ts_mol_name_inact)
            yield i, ts_mol_name_act, ts_mol_name_inact


def get_centroids(cs, df, num):
    return tuple(list(df[df.index == x[0]].values[0]) for x in cs if len(x) >= num)


def trainingset_formation(input_mols, path_ts, mode_train_set, fcfp4, clust_stat, threshold):

    os.makedirs(path_ts, exist_ok=True)
    clust_size, max_num_acts = 5, 5

    if (1 not in mode_train_set) and (2 not in mode_train_set):
        return 'Wrong value of parameter mode_train_set. That should be 1 and/or 2.'

    df_mols, fp = read_file(input_mols, fcfp4)

    list_ts = []
    if 2 in mode_train_set:
        cs = gen_cluster_subset_butina(
            fp,
            threshold
        )
        cs_inact = gen_cluster_subset_butina(
            fp[min(df_mols[df_mols['activity'] == 0].index):],
            threshold
        )
        # get_centroids() returns tuple of tuples with mol names and their SMILES
        centroids_inact = get_centroids(cs_inact, df_mols, clust_size)

        if clust_stat:
            save_cluster_stat(
                cs,
                df_mols,
                df_mols[df_mols['activity'] == 1].index.tolist(),
                clust_stat
            )

        for i, act_ts, inact_ts in diff_binding_mode(
                                        cs, df_mols,
                                        df_mols[df_mols['activity'] == 1].index.tolist(),
                                        centroids_inact, max_num_acts):

            output = os.path.join(path_ts, f't{i}.smi')
            list_ts.append(output)
            save_train_set(output, act_ts, inact_ts)

    if 1 in mode_train_set:
        # process actives
        cs_act = gen_cluster_subset_butina(
            fp[:min(df_mols[df_mols['activity'] == 0].index)],
            threshold
        )
        centroids_act = get_centroids(cs_act, df_mols, clust_size)

        # if number active centroids is less than the minimum number of molecules in the centroid training set
        if len(centroids_act) < max_num_acts:
            return list_ts

        # process inactives
        cs_inact = gen_cluster_subset_butina(
            fp[min(df_mols[df_mols['activity'] == 0].index):],
            threshold
        )
        centroids_inact = get_centroids(cs_inact, df_mols, clust_size)

        output = os.path.join(path_ts, 'centroids.smi')
        list_ts.append(output)
        save_train_set(output, centroids_act, centroids_inact)

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
