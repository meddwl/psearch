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
    """Build the CLI argument parser for select_training_set."""
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
                        help='Training-set strategy: 1 = one set from cluster centroids; '
                             '2 = one set per cluster; supply both (e.g. -ts 1 2) to run both strategies.')
    parser.add_argument('--fcfp4', action='store_true', default=False,
                        help='If set FCFP4 fingerprints will be used for compound clustering, '
                             'otherwise pharmacophore fingerprints will be used.')
    parser.add_argument('-t', '--threshold_clust', metavar='NUMERIC', type=float, default=0.4,
                        help='threshold for clustering data by Butina algorithm')
    parser.add_argument('-s', '--save_statistics', metavar='FILENAME', default=None,
                        help='If a file path is provided, cluster statistics (size, active/inactive counts) '
                             'will be written to that file.')
    return parser


def read_file(fname, fcfp4):
    """Read a tab-separated SMILES file and compute fingerprints for each compound.

    Sorts compounds by activity (actives first) and writes a sorted copy alongside the
    input file (_sorted.smi). Computes either FCFP4 (Morgan radius-2 feature-based) or
    2D pharmacophore fingerprints for use in Butina clustering.

    Args:
        fname: Path to the tab-separated SMILES file (columns: SMILES, mol_id, activity).
        fcfp4: If True, compute FCFP4 fingerprints; otherwise compute 2D pharmacophore fps.

    Returns:
        Tuple of (df, fp_list) where df is the sorted DataFrame and fp_list contains
        one fingerprint object per compound in the same order.
    """
    df = pd.read_csv(fname, sep='\t', header=None, names=['smiles', 'mol_name', 'activity'])
    if not Chem.MolFromSmiles(df.at[0, 'smiles']):
        df = df.drop(index=0)

    df['activity'] = df['activity'].astype(np.int64)
    df = df.sort_values(by=['activity'], ascending=False).reset_index(drop=True)
    df.to_csv(os.path.splitext(fname)[0] + '_sorted.smi', sep='\t', index=None)

    if fcfp4:
        fp = []
        for smiles in df['smiles']:
            mol = Chem.MolFromSmiles(smiles)  # updated: None guard - MolFromSmiles returns None silently on parse failure
            if mol is None:
                raise ValueError(f"Could not parse SMILES: {smiles!r}")
            fp.append(AllChem.GetMorganFingerprint(mol, 2, useFeatures=True))
    else:
        featfactory = load_factory()
        sigfactory = SigFactory(featfactory, minPointCount=2, maxPointCount=3, trianglePruneBins=False)
        sigfactory.SetBins([(0, 2), (2, 5), (5, 8)])
        sigfactory.Init()
        fp = []
        for smiles in df['smiles']:
            mol = Chem.MolFromSmiles(smiles)  # updated: None guard - MolFromSmiles returns None silently on parse failure
            if mol is None:
                raise ValueError(f"Could not parse SMILES: {smiles!r}")
            fp.append(Generate.Gen2DFingerprint(mol, sigfactory))
    return df, fp


def gen_cluster_subset_butina(fps, cutoff):
    """Cluster fingerprints by Tanimoto similarity using the Butina algorithm.

    Args:
        fps: List of fingerprint objects (RDKit Morgan or pharmacophore fps).
        cutoff: Tanimoto distance cutoff (1 − similarity) for cluster membership.

    Returns:
        Tuple of tuples; each inner tuple contains the indices (into fps) of compounds
        in one cluster, with the centroid as the first element.
    """
    dists = []
    for i in range(len(fps)-1):
        sims = DataStructs.BulkTanimotoSimilarity(fps[i], fps[i+1:])
        dists.extend([1 - x for x in sims])
    cs = Butina.ClusterData(dists, len(fps), cutoff, isDistData=True)
    return cs


def save_cluster_stat(cs_index, df, index_acts, clust_stat):
    """Write per-cluster statistics to a file object.

    Each cluster entry records the cluster index, total size, active fraction, and
    the positional indices of all member compounds.

    Args:
        cs_index: Tuple of tuples of compound indices (output of gen_cluster_subset_butina).
        df: DataFrame with compound data (used for index lookup).
        index_acts: List of DataFrame indices corresponding to active compounds.
        clust_stat: Writable file object where statistics will be written.
    """
    clust_stat.write("""#The saved file contains statistical information about the cluster 
#and the index of the molecule according to its location in the input file \n""")
    for i, cluster in enumerate(cs_index):
        i_act = len(set(cluster).intersection(index_acts))
        clust_stat.write(f'\ncluster {i}, cluster length {len(cluster)}, share of active {round(i_act/len(cluster), 3)}\n')
        clust_stat.write(','.join(map(str, cluster)) + '\n')


def save_train_set(output, act_ts, inact_ts):
    """Write active and inactive training-set compounds to a tab-separated .smi file.

    Args:
        output: Path to the output file.
        act_ts: Iterable of (SMILES, mol_id, activity) tuples for actives.
        inact_ts: Iterable of (SMILES, mol_id, activity) tuples for inactives.
    """
    with open(output, 'wt') as f:
        f.write('\n'.join(f'{smiles}\t{mol_name}\t{activity}' for smiles, mol_name, activity in act_ts) + '\n')
        f.write('\n'.join(f'{smiles}\t{mol_name}\t{activity}' for smiles, mol_name, activity in inact_ts) + '\n')


def diff_binding_mode(cs, df, index_acts, inact_centroids, min_num):
    """Yield per-cluster training sets for clusters with sufficient active representation.

    Selects clusters containing at least `min_num` actives and yields active and
    inactive training-set entries for each such cluster. Inactive sets are seeded
    with `inact_centroids` to ensure diverse inactive coverage.

    Args:
        cs: Tuple of tuples of compound indices (output of gen_cluster_subset_butina).
        df: DataFrame with columns smiles, mol_name, activity.
        index_acts: List of DataFrame indices corresponding to active compounds.
        inact_centroids: Array of (SMILES, mol_id, activity) rows for inactive centroids.
        min_num: Minimum number of actives required for a cluster to be included.

    Yields:
        Tuples of (cluster_index, act_ts_array, inact_ts_set).
    """
    for i, c in enumerate(cs):
        if len(set(c).intersection(index_acts)) >= min_num:
            dfc = df.iloc[list(c)]
            ts_mol_name_act = dfc[dfc['activity'] == 1][:5].values
            ts_mol_name_inact = np.append(dfc[dfc['activity'] == 0][:5].values, inact_centroids, axis=0)
            ts_mol_name_inact = set(tuple(element) for element in ts_mol_name_inact)
            yield i, ts_mol_name_act, ts_mol_name_inact


def get_centroids(cs, df, num):
    """Return the centroid compound row from each cluster that contains at least `num` members.

    The centroid is the first element of each cluster tuple (as returned by Butina).

    Args:
        cs: Tuple of tuples of compound indices (output of gen_cluster_subset_butina).
        df: DataFrame with compound data.
        num: Minimum cluster size required to include its centroid.

    Returns:
        Tuple of DataFrame row arrays ([SMILES, mol_id, activity]) for qualifying centroids.
    """
    return tuple(list(df[df.index == x[0]].values[0]) for x in cs if len(x) >= num)


def trainingset_formation(input_mols, path_ts, mode_train_set, fcfp4, clust_stat, threshold):
    """Cluster input molecules and write training set files for pharmacophore model generation.

    Applies Butina clustering to the active and inactive compound sets. Depending on
    `mode_train_set`, produces either a single centroid-based training set (strategy 1),
    per-cluster training sets (strategy 2), or both.

    Args:
        input_mols: Path to the tab-separated SMILES file (columns: SMILES, mol_id, activity).
        path_ts: Directory where training set .smi files will be written.
        mode_train_set: List containing 1 and/or 2 to select the training-set strategy.
        fcfp4: If True, use FCFP4 fingerprints for clustering; otherwise use 2D pharmacophore fps.
        clust_stat: Writable file object for cluster statistics, or None to skip.
        threshold: Butina distance cutoff (Tanimoto dissimilarity) for clustering.

    Returns:
        List of paths to the written training set .smi files, or an error string if
        mode_train_set contains neither 1 nor 2.
    """
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
    """CLI entry point for select_training_set: parse arguments and call trainingset_formation."""
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
