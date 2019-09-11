#!/usr/bin/env python3
# author          : Alina Kutlushina
# date            : 01.05.2018
# license         : BSD-3
#==============================================================================

import sys
import argparse
import pandas as pd


def main(in_fname, out_act_fname, out_inact_fname, act_threshold, inact_threshold, label):

    if label:
        df = pd.read_csv(in_fname, sep='\t', header=None)
        df_act = df[df[2] == 'active']
        df_act.to_csv(out_act_fname, sep='\t', index=None, header=None)
        df_inact = df[df[2] == 'inactive']
        df_inact.to_csv(out_inact_fname, sep='\t', index=None, header=None)

        n_act = df_act.shape[0]
        n_inact = df_inact.shape[0]

    else:
        fact = open(out_act_fname, 'wt')
        finact = open(out_inact_fname, 'wt')

        n_act = 0
        n_inact = 0

        try:
            with open(in_fname) as f:
                for line in f:
                    row = line.strip().split(',')
                    if float(row[2]) >= act_threshold:
                        fact.write('\t'.join(row[:3]) + '\n')
                        n_act += 1
                    elif float(row[2]) <= inact_threshold:
                        finact.write('\t'.join(row[:3]) + '\n')
                        n_inact += 1

        finally:
            fact.close()
            finact.close()

    sys.stderr.write('actives: %i, inactives: %i.\n' % (n_act, n_inact))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='select active and inactive compounds'
                                                 'based on given values (act_threshold and inact_threshold)')
    parser.add_argument('-i', '--in', metavar='input.smi', required=True,
                        help='input SMILES file name. It should contain three columns separated by whitespaces: '
                             'SMILES, name, activity. No header.')
    parser.add_argument('-oa', '--out_act', metavar='active.smi', required=True,
                        help='output SMILES file name for active compounds.')
    parser.add_argument('-oi', '--out_inact', metavar='inactive.smi', required=True,
                        help='output SMILES file name for inactive compounds.')
    parser.add_argument('-ta', '--act_threshold', metavar='VALUE', default=8.0,
                        help='specify threshold used to determine active compounds.'
                             'Compounds having activity higher or equal to the given'
                             'value will be recognized as active. Default: 8.')
    parser.add_argument('-ti', '--inact_threshold', metavar='inact_threshold', default=6.0,
                        help='specify threshold used to determine inactive compounds.'
                             'Compounds having activity less or equal to the given'
                             'value will be recognized as inactive. Default: 6.')
    parser.add_argument('-l', '--label', action='store_true', default=False,
                        help='a criterion of molecules separation. '
                             'If True - absolute separation of molecules into active and inactive.'
                             'If False - separation of molecules into active and inactive by value')

    args = vars(parser.parse_args())
    for o, v in args.items():
        if o == "in": in_fname = v
        if o == "out_act": out_act_fname = v
        if o == "out_inact": out_inact_fname = v
        if o == "act_threshold": act_threshold = float(v)
        if o == "inact_threshold": inact_threshold = float(v)
        if o == "label": label = v

    main(in_fname=in_fname,
         out_act_fname=out_act_fname,
         out_inact_fname=out_inact_fname,
         act_threshold=act_threshold,
         inact_threshold=inact_threshold,
         label=label)
