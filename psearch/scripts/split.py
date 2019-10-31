#!/usr/bin/env python3
# author          : Alina Kutlushina
# date            : 01.05.2018
# license         : BSD-3
#==============================================================================

import sys
import argparse
import pandas as pd


def main(in_fname, out_act_fname, out_inact_fname):
    """
    split a dataset into an active and an inactive sets by status column
    :param in_fname: input .smi file
    :param out_act_fname: path where wlii saved an active set
    :param out_inact_fname: path where will saved an inactive set
    :return: None
    """

    df = pd.read_csv(in_fname, sep=',')

    # df = df[['standardized_canonical_smiles', 'cmp', 'status']]
    df_act = df[df['status'] == 'active']
    df_act.to_csv(out_act_fname, sep='\t', index=None, header=None)
    df_inact = df[df['status'] == 'inactive']
    df_inact.to_csv(out_inact_fname, sep='\t', index=None, header=None)

    sys.stderr.write('actives: %i, inactives: %i.\n' % (df_act.shape[0], df_inact.shape[0]))


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


    args = vars(parser.parse_args())
    for o, v in args.items():
        if o == "in": in_fname = v
        if o == "out_act": out_act_fname = v
        if o == "out_inact": out_inact_fname = v

    main(in_fname=in_fname,
         out_act_fname=out_act_fname,
         out_inact_fname=out_inact_fname)
