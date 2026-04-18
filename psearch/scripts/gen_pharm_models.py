#!/usr/bin/env python3
# author          : Alina Kutlushina
# date            : 10.01.2019
# license         : BSD-3
# ==============================================================================

import os.path
import sys
import time
import argparse
import pandas as pd
# from psearch.database import DB
from psearch.database import load_model
from pmapper.pharmacophore import Pharmacophore


def _keep_best_models(df, df_sub, save_files, nfeatures):
    """Filter df_sub to retain only pharmacophore models whose hashes appear in df (the high-precision set).

    Optionally writes internal statistics and the filtered model table to disk.

    Args:
        df: DataFrame of high-precision model hashes (output of calc_internal_stat).
        df_sub: Full candidate model DataFrame to filter.
        save_files: Either False (no saving) or a list [output_dir, cluster_id] for file output.
        nfeatures: Current feature count, used to label output filenames.

    Returns:
        Filtered df_sub containing only models present in df.
    """
    df_sub = df_sub[df_sub['hash'].isin(set(df['hash']))].reset_index(drop=True)
    if save_files:
        df.to_csv(os.path.join(save_files[0], f'internal_statistics-{save_files[1]}-f{nfeatures}.txt'), index=None, sep='\t')
        df_sub.to_csv(os.path.join(save_files[0], f'pharmacophore_models-{save_files[1]}-f{nfeatures}.txt'), index=None, sep='\t')
    return df_sub


def _gen_quadruplets(db_path, pp_train_set, lower, tol, bin_step):
    """Enumerate all `lower`-feature pharmacophore sub-graphs for each conformer in the training set.

    Iterates over all stereoisomers and conformers of each training-set compound, builds
    Pharmacophore objects, and yields one record per sub-graph enumerated by pmapper.

    Args:
        db_path: Path to the psearch SQLite database file.
        pp_train_set: Path to the training set .smi file (columns: SMILES, mol_id, activity).
        lower: Number of features to enumerate (typically 4 for quadruplets).
        tol: Stereocentre sign tolerance passed to pmapper's iterate_pharm.
        bin_step: Bin width (Å) for Pharmacophore objects.

    Yields:
        Tuples of (activity, mol_name, isomer_id, conf_id, hash, labels).
    """
    train_set_list = [name.strip().split() for name in open(pp_train_set).readlines()]
    for _, mol_name, activity in train_set_list:
        try:
            # dict_coords = db.get_pharm(mol_name)
            dict_coords = load_model(db_path, mol_name).pharm_dict
        except KeyError:
            raise KeyError(f"Molecule '{mol_name}' not found in database")
        for isomer_id, list_coords in dict_coords.items():
            for conf_id, coord in enumerate(list_coords):
                pharm = Pharmacophore(bin_step=bin_step, cached=True)
                pharm.load_from_feature_coords(coord)
                if pharm:
                    for hash, labels in pharm.iterate_pharm(lower, lower, tol):
                        yield activity, mol_name, isomer_id, conf_id, hash, labels


def _plus_one_feature(db_path, df_sub, bin_step):
    """Extend existing pharmacophore models by one additional feature.

    For each (mol, stereo, conformer) group in df_sub, retrieves the full pharmacophore
    from the database and uses pmapper's iterate_pharm1 to enumerate all extensions of
    the current feature sets by exactly one additional feature.

    Args:
        db_path: Path to the psearch SQLite database file.
        df_sub: DataFrame of current-complexity candidate models (from _keep_best_models).
        bin_step: Bin width (Å) for Pharmacophore objects.

    Yields:
        Tuples of (activity, mol_name, isomer_id, conf_id, hash, labels).
    """
    df_sub_group = df_sub.groupby(['activity', 'mol_name', 'isomer_id', 'conf_id'])
    cols = ['activity', 'mol_name', 'isomer_id', 'conf_id']
    pharm_cache = {}  # avoid repeated DB lookups for the same conformer
    for df_group in df_sub_group:
        activity, mol_name, isomer_id, conf_id = df_group[1][cols].iloc[0]
        label_ids = [tuple(map(int, lbls.split(','))) for lbls in df_group[1]['feature_ids'].tolist()]
        cache_key = (mol_name, isomer_id, conf_id)
        if cache_key not in pharm_cache:
            # pharm_cache[cache_key] = db.get_pharm(mol_name)[isomer_id][conf_id]
            pharm_cache[cache_key] = load_model(db_path, mol_name).pharm_dict[isomer_id][conf_id]
        pharm = Pharmacophore(bin_step=bin_step, cached=True)
        pharm.load_from_feature_coords(pharm_cache[cache_key])
        if pharm:
            for hash, labels in pharm.iterate_pharm1(label_ids):
                yield activity, mol_name, isomer_id, conf_id, hash, labels


def gen_models(def_generator):
    """Collect pharmacophore sub-graph records from a generator and return a summary DataFrame.

    Aggregates per-compound hit counts for each hash (pharmacophore) and returns a
    DataFrame used by calc_internal_stat to evaluate model quality.

    Args:
        def_generator: Generator yielding (activity, mol_name, isomer_id, conf_id, hash, labels)
                       tuples — typically _gen_quadruplets or _plus_one_feature.

    Returns:
        DataFrame with columns [activity, hash, count, mol_name, isomer_id, conf_id, feature_ids],
        sorted by activity and descending hit count. Empty DataFrame if no records.
    """
    data = []
    for activity, mol_name, isomer_id, conf_id, hash, labels in def_generator:
        data.append([activity, mol_name, isomer_id, conf_id, hash, ','.join(map(str, labels))])
    df = pd.DataFrame(data, columns=['activity', 'mol_name', 'isomer_id', 'conf_id', 'hash', 'feature_ids'])
    if df.empty:
        return df
    count_df = df.drop_duplicates(subset=['activity', 'mol_name', 'hash'])
    count_df = count_df.groupby(['activity', 'hash'], sort=True).size().reset_index(name='count')
    df = pd.merge(df, count_df, on=['activity', 'hash'], how='outer')
    df = df.sort_values(by=['activity', 'count', 'hash'], ascending=False)
    return df[['activity', 'hash', 'count', 'mol_name', 'isomer_id', 'conf_id', 'feature_ids']]


def strategy_extract_trainset(df, clust_strategy):
    """Select the top-performing pharmacophore models according to the clustering strategy.

    Strategy 1 (centroid training set): keeps models with F0.5 >= 0.8 (precision-biased).
    Strategy 2 (per-cluster training sets): keeps models with F2 >= 0.8 (recall-biased),
    favouring coverage of actives across diverse binding modes.

    Args:
        df: DataFrame with columns [hash, TP, FP, precision, recall, F2, F05].
        clust_strategy: 1 for centroid-based selection; 2 for per-cluster selection.

    Returns:
        Filtered DataFrame containing only the selected high-quality models.
    """
    if clust_strategy == 2:
        df = df.sort_values(by=['recall', 'F2', 'F05'], ascending=False).reset_index(drop=True)
        if df['F2'].iloc[0] == 1.0:
            df = df[(df['recall'] == 1.0) & (df['F2'] == 1.0)]
        elif df[df['F2'] >= 0.8].shape[0] <= 100:
            df = df[(df['recall'] == 1) & (df['F2'] >= 0.8)]
        else:
            df = df[(df['recall'] == 1) & (df['F2'] >= df['F2'].loc[100])]
    elif clust_strategy == 1:
        df = df.sort_values(by=['recall', 'F05', 'F2'], ascending=False).reset_index(drop=True)
        df = df[df['F05'] >= 0.8] if df[df['F05'] >= 0.8].shape[0] <= 100 else df[df['F05'] >= df['F05'].loc[100]]
    return df


def calc_internal_stat(df, positives, clust_strategy, designating):
    """Compute internal validation statistics for candidate pharmacophore models.

    Calculates TP, FP, precision, recall, F2, and F0.5 for each unique pharmacophore
    hash and then calls strategy_extract_trainset to retain only top-quality models.

    Args:
        df: Deduplicated DataFrame with columns [activity, hash, count].
        positives: Total number of active compounds in the training set.
        clust_strategy: Passed to strategy_extract_trainset (1 = centroid, 2 = per-cluster).
        designating: List of two activity label strings, e.g. ['1', '0'] for active/inactive.

    Returns:
        DataFrame with columns [hash, TP, FP, precision, recall, F2, F05] for selected models.
    """
    if df[df['activity'] == designating[1]].empty:
        df['FP'] = [0] * df.shape[0]
        df = df.rename(columns={'count': 'TP'})
    else:
        df = df[df['activity'] == designating[0]].rename(columns={'count': 'TP'}).merge(
             df[df['activity'] == designating[1]].rename(columns={'count': 'FP'}),
             on='hash', how='outer')
        df.loc[df['FP'].isnull(), 'FP'] = 0
        df.loc[df['TP'].isnull(), 'TP'] = 0
    df = df[['hash', 'TP', 'FP']]
    df['precision'] = round(df['TP'] / (df['TP'] + df['FP']), 3)
    df['recall'] = round(df['TP'] / positives, 3)
    df['F2'] = round(5 * ((df['precision'] * df['recall']) / (4 * df['precision'] + df['recall'])), 3)
    df['F05'] = round(1.25 * ((df['precision'] * df['recall']) / (0.25 * df['precision'] + df['recall'])), 3)
    df = df[['hash', 'TP', 'FP', 'precision', 'recall', 'F2', 'F05']]
    # difference ways to check out the best models
    df = strategy_extract_trainset(df, clust_strategy)
    return df


def save_models_xyz(db_path, db_name, df_sub, path_pma, bin_step, cluster_id, num_ids):
    """Write pharmacophore models to .xyz files in path_pma.

    Each unique model hash in df_sub is written as a separate .xyz file using the
    representative conformer stored in the database. Files are named
    {db_name}.{cluster_id}_f{num_ids}_p{index}.xyz.

    Args:
        db_path: Path to the psearch SQLite database file.
        db_name: Database name string used as a filename prefix.
        df_sub: DataFrame of selected models (output of _keep_best_models).
        path_pma: Directory where .xyz files will be written.
        bin_step: Bin width (Å) for Pharmacophore objects.
        cluster_id: Cluster identifier string used in filenames.
        num_ids: Feature count used in filenames.

    Returns:
        Number of model files written.
    """
    data = df_sub.drop_duplicates(subset=['hash']).values
    for num, (_, hash, count, mol_name, isomer_id, conf_id, feature_ids) in enumerate(data):
        pharm = Pharmacophore(bin_step=bin_step, cached=True)
        # pharm.load_from_feature_coords(db.get_pharm(mol_name)[isomer_id][conf_id])
        pharm.load_from_feature_coords(load_model(db_path, mol_name).pharm_dict[isomer_id][conf_id])
        pharm.save_to_xyz(os.path.join(path_pma, f"{db_name}.{cluster_id}_f{num_ids}_p{num}.xyz"),
                          tuple(map(int, feature_ids.split(','))))
    return len(data)


def gen_pharm_models(in_db, out_pma, trainset, tolerance, bin_step, current_nfeatures, upper, nfeatures, save_statistics):
    """Iteratively generate ligand-based pharmacophore models from a training set.

    Starts from `current_nfeatures`-feature pharmacophore seeds (quadruplets or higher),
    computes internal statistics, prunes low-quality candidates, then extends survivors
    by one feature at a time until `upper` features or no further improvement.

    Args:
        in_db: Path to the psearch database (.dat file).
        out_pma: Directory where final .xyz model files will be written.
        trainset: Path to the training set .smi file (columns: SMILES, mol_id, activity).
        tolerance: Stereocentre sign tolerance for pmapper feature enumeration.
        bin_step: Bin width (Å) for pharmacophore discretisation.
        current_nfeatures: Starting feature count for model generation (typically 4).
        upper: Maximum feature count; None for maximum possible complexity.
        nfeatures: Minimum feature count at which models are saved; None saves only the final set.
        save_statistics: If True, write intermediate statistics to an 'intermediate_data' subdirectory.
    """
    time_start = time.time()
    db_name = os.path.splitext(os.path.basename(in_db))[0]
    os.makedirs(out_pma, exist_ok=True)
    cluster_id = os.path.splitext(os.path.basename(trainset))[0]
    designating = ['1', '0']  # molecular activity
    clust_strategy = 1 if cluster_id == 'centroids' else 2
    positives = len([line for line in open(trainset).readlines() if line.strip().split()[2] == designating[0]])
    # db = DB(in_db, flag='r')
    df_sub = gen_models(_gen_quadruplets(in_db, trainset, current_nfeatures, tolerance, bin_step))
    df = calc_internal_stat(df_sub[['activity', 'hash', 'count']].drop_duplicates(subset=['activity', 'hash']),
                            positives, clust_strategy, designating)
    if df.empty:
        sys.stderr.write(f'train set {cluster_id}: no {current_nfeatures}-points pharmacophore models\n')
        return

    if save_statistics:
        path_files = os.path.join(out_pma, 'intermediate_data')
        os.makedirs(path_files, exist_ok=True)
        save_statistics = [path_files, cluster_id]

    df_sub = _keep_best_models(df, df_sub, save_statistics, current_nfeatures)
    if nfeatures is not None:
        if current_nfeatures >= nfeatures:
            _ = save_models_xyz(in_db, db_name, df_sub[df_sub['activity'] == designating[0]],
                                out_pma, bin_step, cluster_id, current_nfeatures)

    while True:
        if current_nfeatures == upper:
            break
        current_nfeatures += 1
        df_sub_2 = gen_models(_plus_one_feature(in_db, df_sub, bin_step))
        df = calc_internal_stat(df_sub_2[['activity', 'hash', 'count']].drop_duplicates(subset=['activity', 'hash']),
                            positives, clust_strategy, designating)
        if df.empty:
            break

        df_sub = _keep_best_models(df, df_sub_2, save_statistics, current_nfeatures)
        if nfeatures is not None:
            if current_nfeatures >= nfeatures:
                _ = save_models_xyz(in_db, db_name, df_sub[df_sub['activity'] == designating[0]],
                                    out_pma, bin_step, cluster_id, current_nfeatures)

    num_models = save_models_xyz(in_db, db_name, df_sub[df_sub['activity'] == designating[0]], out_pma, bin_step, cluster_id, current_nfeatures)
    sys.stderr.write(f'train set {cluster_id}: {num_models} models ({round(time.time()-time_start, 3)}s)\n')
    sys.stderr.flush()


def create_parser():
    """Build the CLI argument parser for gen_pharm_models."""
    parser = argparse.ArgumentParser(description='Iteratively create ligand-based pharmacophore models.',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--database', metavar='FILENAME.dat', required=True,
                        help='Input psearch database file (.dat) with precomputed conformers and pharmacophores.')
    parser.add_argument('-o', '--models', metavar='DIRNAME', required=False, default=None,
                        help='Output path to a folder where will be saved the created pharmacophore models.'
                             'If omitted, the path will be generated automatically relative to project directory.')
    parser.add_argument('-ts', '--trainset', metavar='FILENAME.txt', required=True,
                        help='Path to tab-separated txt file with information about molecules from a training set.'
                             'Required columns: SMILES, MOL_ID, ACTIVITY')
    parser.add_argument('-tol', '--tolerance', metavar='NUMERIC', type=float, default=0,
                        help='tolerance used for calculation of a stereoconfiguration sign.')
    parser.add_argument('-b', '--bin_step', metavar='INTEGER', type=int, default=1,
                        help='Bin width (Å) for discretising pharmacophore feature coordinates. '
                             'Must match the value used when building the database with gen_db.')
    parser.add_argument('-u', '--upper', type=int,  metavar='INTEGER', default=None,
                        help='Maximum number of pharmacophore features per model. '
                             'If omitted, models up to the maximum possible complexity are generated.')
    parser.add_argument('-l', '--lower', metavar='INTEGER', type=int, default=3,
                        help='Minimum number of pharmacophore features to include when generating models. '
                             'Models with fewer features are not created.')
    parser.add_argument('-f', '--save_model_complexity', type=int, metavar='INTEGER', default=None,
                        help='Minimum feature count at which models are written to disk. '
                             'If omitted, only the most complex models are saved.')
    parser.add_argument('-s', '--save_statistics', action='store_true', default=False,
                        help='Save intermediate per-feature-count model statistics to an intermediate_data subdirectory. '
                             'Useful for diagnosing model generation quality at each step.')
    return parser


if __name__ == '__main__':
    parser = create_parser()
    args = parser.parse_args()
    gen_pharm_models(in_db=os.path.abspath(args.database),
                     trainset=os.path.abspath(args.trainset),
                     out_pma=os.path.abspath(args.models) if args.models else os.path.join(os.path.split(os.path.abspath(args.database))[0], 'models'),
                     bin_step=int(args.bin_step),
                     tolerance=args.tolerance,
                     current_nfeatures=args.lower,
                     upper=int(args.upper) if args.upper is not None else None,
                     nfeatures=int(args.save_model_complexity) if args.save_model_complexity is not None else None,
                     save_statistics=args.save_statistics)
