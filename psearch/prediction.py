#!/usr/bin/env python3
# author          : Alina Kutlushina
# date            : 02.12.20
# license         : BSD-3
# ==============================================================================

__author__ = 'Alina Kutlushina'

import os
import pandas as pd
import argparse
from collections import defaultdict
default_modelstat = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'pharmacophores', 'pharmacophores_stat.csv')

def calc_probability(df_vs, df_precision, target_id, scoring_scheme):
    """Compute the predicted activity probability of hits agains a target from VS results.

    Multiplies each model's hit vector by its precision, then collapses across models
    using the chosen scoring scheme.

    Args:
        df_vs: DataFrame (index=model_id, columns=mol_id) with 1/NaN hit indicators.
        df_precision: DataFrame (index=model_id) with a 'precision' column.
        target_id: Target identifier used to label the result row.
        scoring_scheme: 'mean' averages probabilities across models;
                        'max' takes the highest individual model probability.

    Returns:
        Series of rounded (3 d.p.) activity probabilities indexed by mol_id.
    """
    df = df_vs.mul(df_precision['precision'].astype(float), axis=0)
    if scoring_scheme == 'mean':
        df.loc[target_id] = df.mean(axis=0, skipna=True)
    elif scoring_scheme == 'max':
        df.loc[target_id] = df.max(axis=0, skipna=True)
    return round(df.loc[target_id], 3)


def input_processing(list_vs, models_list):
    """Build a hit-indicator DataFrame from VS result files.

    Reads each hit-list file in list_vs, extracts the model_id from the filename
    (format: target_id.model_id.txt), and marks each hit compound with 1. Compounds
    not hit by a model appear as NaN.

    Args:
        list_vs: List of paths to screen_db hit-list .txt files for one target.
        models_list: List of model_ids to use as the DataFrame index.

    Returns:
        DataFrame (index=model_id, columns=mol_id) with 1/NaN hit indicators.
    """
    data = {}
    for ff in list_vs:
        ph = os.path.splitext(os.path.basename(ff))[0].split('.')[1]
        mols = [i.strip().split()[0] for i in open(ff).readlines()]
        for mol_id in mols:
            if mol_id not in data:
                data[mol_id] = {}
            data[mol_id][ph] = 1
    df = pd.DataFrame(data, index=models_list)
    return df


def main(pp_vs, pp_models_stat, scoring_scheme, pp_output):
    """Load VS results and model statistics, compute per-target probabilities, and save output.

    Args:
        pp_vs: Path to the directory containing screen_db hit-list .txt files.
        pp_models_stat: Path to a TSV file with columns target_id, model_id, precision.
        scoring_scheme: 'mean' or 'max' — see calc_probability.
        pp_output: Path to the output TSV file where predictions will be written.
    """
    df_models_stat = pd.read_csv(pp_models_stat, sep='\t', index_col='model_id')
    vs_files = defaultdict(list)
    for fname in os.listdir(pp_vs):
        if not fname.endswith('.txt'):   # skip .sdf files and other non-result files
            continue
        target_id = os.path.splitext(fname)[0].split('.')[0]
        vs_files[target_id].append(os.path.join(pp_vs, fname))

    result_series = []
    target_ids = sorted(vs_files.keys())
    for target_id in target_ids:
        df_models = df_models_stat[df_models_stat['target_id'] == target_id]
        df_vs = input_processing(vs_files[target_id], df_models.index.tolist())
        df_res = calc_probability(df_vs, df_models, target_id, scoring_scheme)
        result_series.append(df_res)

    if not result_series:
        pd.DataFrame().to_csv(pp_output, sep='\t')
        return

    res = pd.concat(result_series, axis=1)
    res.index.name = 'mol_id'
    res = res.sort_values(by=res.columns.tolist()[0], ascending=False)
    res.to_csv(pp_output, sep='\t')


def entry_point():
    """CLI entry point for prediction: parse arguments and call main."""
    parser = argparse.ArgumentParser(description='Determination of the probability of activity of a molecule(-s)'
                                                 'based on pharmacophore VS result',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-s', '--path_vs', metavar='DIRNAME', required=True,
                        help='Path to the directory containing screen_db hit-list .txt files '
                             '(one file per pharmacophore model, named target_id.model_id.txt).')
    parser.add_argument('-p', '--pharm_stat', metavar='FILENAME', default=None,
                        help='file with the calculated precision of pharmacophore models. '
                             'By default, statistics of psearch pharmacophore models are used.'
                             'Required headers: "target_id", "model_id", "precision"')
    parser.add_argument('-f', '--scoring_scheme', metavar='KEYWORD', default='mean',
                        help="Consensus scoring scheme: 'max' uses the highest individual model probability; "
                             "'mean' averages probabilities across all models for a target.")
    parser.add_argument('-o', '--output', metavar='FILENAME', default=None,
                        help='Output tab-separated text file where per-target activity predictions will be written.')

    args = parser.parse_args()
    if args.scoring_scheme not in ('mean', 'max'):
        parser.error(f"--scoring_scheme must be 'mean' or 'max', got '{args.scoring_scheme}'")
    output = os.path.abspath(args.output) if args.output else os.path.join(os.path.abspath(args.path_vs), 'prediction.txt')
    os.makedirs(os.path.dirname(output), exist_ok=True)
    pharm_stat = os.path.abspath(args.pharm_stat) if args.pharm_stat else default_modelstat
    main(os.path.abspath(args.path_vs), pharm_stat, args.scoring_scheme, output)


if __name__ == '__main__':
    entry_point()