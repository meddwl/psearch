#!/usr/bin/env python3
# author          : Alina Kutlushina
# date            : 01.05.2018
# license         : BSD-3
# ==============================================================================

import os
import sys
import time
import argparse
import pandas as pd
from rdkit import Chem
from pmapper.pharmacophore import Pharmacophore as P


def max_edge(model):
    """Return the length (Å) of the longest pairwise distance between any two features in a pharmacophore model.

    Args:
        model: Path to a .xyz pharmacophore model file.

    Returns:
        Float — maximum inter-feature Euclidean distance in Å.
    """
    p = P()
    p.load_from_xyz(model)
    coords = p.get_feature_coords()
    edge = 0
    for i, c1 in enumerate(coords):
        for j, c2 in enumerate(coords[i + 1:]):
            e = ((c1[1][0] - c2[1][0]) ** 2 + (c1[1][1] - c2[1][1]) ** 2 + (c1[1][2] - c2[1][2]) ** 2) ** (1 / 2)
            if e > edge:
                edge = e
    return edge


def get_external_stat(path_mols, path_ts, path_pma, pp_screen, model_id):
    """Compute external validation statistics for a single pharmacophore model.

    Compares virtual screening hits against the external test set (compounds excluded from
    the training set) to compute TP, FP, precision, recall, F-scores, balanced accuracy,
    and enrichment factor.

    Args:
        path_mols: Path to the full tab-separated SMILES file (columns: SMILES, mol_id, activity).
        path_ts: Path to the training set .smi file for this model's cluster.
        path_pma: Path to the .xyz pharmacophore model file.
        pp_screen: Path to the screen_db hit-list .txt file for this model.
        model_id: Model identifier string (used in the returned results tuple).

    Returns:
        Tuple of (target_id, model_id, TP, FP, P, N, precision, recall, FPR,
                  F1, F2, F05, BA, EF, uniq_features, max_dist, feature_labels).
    """
    target_id = os.path.splitext(os.path.basename(path_mols))[0]
    p = P()
    p.load_from_xyz(path_pma)
    feature_coords = p.get_feature_coords()
    labels = ''.join(label for label, xyz in feature_coords)
    num_uniq_features = len(set(xyz for label, xyz in feature_coords))
    medge = max_edge(path_pma)

    ts_act_mol = []
    ts_inact_mol = []
    for ii in open(path_ts).readlines():
        line = ii.strip().split()
        if line[-1] == '1':
            ts_act_mol.append(line[1])
        else:
            ts_inact_mol.append(line[1])

    df_mols = pd.read_csv(path_mols, sep='\t', header=None).rename(columns={0: 'smiles', 1: 'mol_name', 2: 'activity'})
    if not Chem.MolFromSmiles(df_mols.at[0, 'smiles']):
        df_mols.drop(index=0, inplace=True)
    df_mols = df_mols.astype(str)
    df_act = df_mols[(df_mols['activity'] == '1') & (~df_mols['mol_name'].isin(ts_act_mol))]
    df_inact = df_mols[(df_mols['activity'] == '0') & (~df_mols['mol_name'].isin(ts_inact_mol))]

    if os.path.exists(pp_screen):
        res_screen = [ii.strip().split()[0] for ii in open(pp_screen).readlines()]
        act_screen = set(res_screen) & set(df_act['mol_name'])
        inact_screen = set(res_screen) & set(df_inact['mol_name'])
    else:
        act_screen = []
        inact_screen = []

    p = df_act.shape[0]
    n = df_inact.shape[0]
    tp = len(act_screen)
    fp = len(inact_screen)
    tn = n - fp

    try:
        recall = tp / p
    except ZeroDivisionError:
        recall = None
    try:
        tnr = tn / n
        fpr = fp / n
    except ZeroDivisionError:
        tnr, fpr = None, None
    if recall and tnr:
        ba = (recall + tnr) / 2
    else:
        ba = None

    try:
        precision = tp / (tp + fp)
    except ZeroDivisionError:
        precision = None

    if precision:
        ef = precision / (p / (p + n))
        f1 = (2 * precision * recall) / (precision + recall)
        f2 = (5 * precision * recall) / (4 * precision + recall)
        f05 = (1.25 * precision * recall) / (0.25 * precision + recall)
    else:
        ef, f1, f2, f05 = None, None, None, None
    return target_id, model_id, tp, fp, p, n, precision, recall, fpr, f1, f2, f05, ba, ef, num_uniq_features, medge, labels


def calc_stat(path_mols, path_ts, pp_models, path_screen, out_external):
    """Compute external validation statistics for all pharmacophore models in a directory.

    Iterates over .pma/.xyz files in `pp_models`, calls get_external_stat for each,
    and writes a sorted TSV summary to `out_external`.

    Args:
        path_mols: Path to the full tab-separated SMILES file (columns: SMILES, mol_id, activity).
        path_ts: Path to the directory containing training set .smi files.
        pp_models: Path to the directory containing pharmacophore model files (.pma or .xyz).
        path_screen: Path to the directory containing screen_db hit-list .txt files.
        out_external: Path to the output TSV file where validation statistics will be written.
    """
    start_time = time.time()
    os.makedirs(os.path.dirname(out_external), exist_ok=True)
    df_result = pd.DataFrame(columns=['target_id', 'model_id', 'TP', 'FP', 'P', 'N', 'precision', 'recall', 'FPR',
                                      'F1', 'F2', 'F05', 'BA', 'EF', 'uniq_features', 'max_dist', 'features'])
    for enum, fmodel in enumerate(sorted(os.listdir(pp_models))):
        if os.path.isfile(os.path.join(pp_models, fmodel)) and (fmodel.endswith('.pma') or fmodel.endswith('.xyz')):
            target_id, model_id = os.path.splitext(fmodel)[0].split('.')
            results = get_external_stat(path_mols=path_mols,
                                        path_ts=os.path.join(path_ts, f'{model_id.split("_")[0]}.smi'),
                                        path_pma=os.path.join(pp_models, fmodel),
                                        pp_screen=os.path.join(path_screen, f'{target_id}.{model_id}.txt'),
                                        model_id=model_id)
            if results:
                df_result.loc[enum] = results
            else:
                continue

    df_result = df_result.sort_values(by=['recall', 'F05', 'F2'], ascending=False)
    df_result = round(df_result, 3)
    df_result.to_csv(out_external, index=None, sep='\t')
    sys.stderr.write(f'\n{os.path.basename(out_external)}: ({round(time.time() - start_time, 3)}s)\n')


def create_parser():
    """Build the CLI argument parser for external_statistics."""
    parser = argparse.ArgumentParser(description='External statistics calculation. '
                                                 'If the metric cannot be calculated (dividing by zero), '
                                                 'NaN will be printed',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--molecules', metavar='FILENAME.smi', required=True,
                        help='The script takes as input a tab-separated SMILES file containing `SMILES`, '
                             '`compound id`, `activity` columns. '
                             'The third column should contain a word 1 or 0. 1 is for actives, 0 is for inactive ones.')
    parser.add_argument('-t', '--trainset', metavar='DIRNAME', required=True,
                        help='Path to the folder containing the pre-computed training set files '
                             '(output of select_training_set).')
    parser.add_argument('-m', '--models', metavar='DIRNAME', required=True,
                        help='Path to the folder containing the pharmacophore model files (.pma or .xyz).')
    parser.add_argument('-s', '--screen', metavar='DIRNAME', required=True,
                        help='path to the folder with the virtual screening results')
    parser.add_argument('-o', '--output', metavar='FILENAME', default=None,
                        help='An output text file where will be saved validation statistics.'
                             'If omitted, the path will be generated automatically relative to project directory.')
    return parser


def entry_point():
    """CLI entry point for external_statistics: parse arguments and call calc_stat."""
    parser = create_parser()
    args = parser.parse_args()
    calc_stat(path_mols=os.path.abspath(args.molecules),
              path_ts=os.path.abspath(args.trainset),
              pp_models=os.path.abspath(args.models),
              path_screen=os.path.abspath(args.screen),
              out_external=os.path.abspath(args.output) if args.output else os.path.join(os.path.dirname(args.molecules), 'result.txt'))


if __name__ == '__main__':
    entry_point()
