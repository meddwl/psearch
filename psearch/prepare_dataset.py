#!/usr/bin/env python3
# author          : Alina Kutlushina
# date            : 01.05.2018
# license         : BSD-3
#==============================================================================
import os
import sys
import time

import pandas as pd
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from multiprocessing import Process

from scripts.gen_db import create_db


def create_parser():
    parser = ArgumentParser(description='', formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--input', metavar='input.smi', type=str, required=True,
                        help='input smi file or multiple files')
    parser.add_argument('-f', '--rdkit_factory', metavar='features.fdef', default=None,
                        help='text file with definition of pharmacophore features in RDKit format. If file name is not '
                             'specified the default file from the script dir will be used. This option has '
                             'a priority over smarts_features.')
    parser.add_argument('-ns', '--nstereo', metavar='stereo_number', default=3,
                        help='number of generated stereoisomers.')
    parser.add_argument('-nc', '--nconf', metavar='conf_number', default=100,
                        help='number of generated conformers.')
    parser.add_argument('-e', '--energy_cutoff', metavar='100', default=100,
                        help='conformers with energy difference from the lowest one greater than the specified '
                             'value will be discarded.')
    parser.add_argument('-r', '--rms', metavar='rms_threshold', default=0.5,
                        help='only conformers with RMS higher then threshold will be kept.')
    parser.add_argument('-b', '--bin_step', default=1,
                        help='binning step. Default: 1.')
    parser.add_argument('-c', '--ncpu', metavar='cpu_number', default=1,
                        help='number of cpus to use for processing of actives and inactives separately. ')
    return parser


def split(in_fname, out_act_fname, out_inact_fname):
    """
    split a dataset into an active and an inactive sets by status column
    :param in_fname: input .smi file
    :param out_act_fname: path where an active set will be saved
    :param out_inact_fname: path where an inactive set will be saved
    :return: None
    """

    df = pd.read_csv(in_fname, sep='\t', header=None)
    df_act = df[df[2] == 'active']
    df_act.to_csv(out_act_fname, sep='\t', index=None, header=None)
    df_inact = df[df[2] == 'inactive']
    df_inact.to_csv(out_inact_fname, sep='\t', index=None, header=None)

    sys.stderr.write('actives: %i, inactives: %i.\n' % (df_act.shape[0], df_inact.shape[0]))


def common(input_fname, db_fname, nstereo, nconf, energy, rms, rdkit_factory, bin_step, ncpu, set_name):
    start = time.time()

    create_db(in_fname=input_fname,
              out_fname=db_fname,
              rdkit_factory=rdkit_factory,
              id_field_name=None,
              nconf=nconf,
              nstereo=nstereo,
              energy=energy,
              rms=rms,
              ncpu=ncpu,
              bin_step=bin_step,
              seed=-1,
              verbose=True)

    sys.stderr.write('prepare {} dataset ({}s)'.format(set_name, time.time() - start))


def main(in_fname, rdkit_factory, nstereo, nconf, energy, rms, bin_step, ncpu):
    """
    launches the entire cycle of data preprocessing: generation of stereoisomers, conformers and a database
    :param in_fname: input .smi file containing information about SMILES, compounds id and its activity status
    :param split_dataset: if True will splited input dasets into active and inactive sets else will not
    :param rdkit_factory: text file with definition of pharmacophore features in RDKit format.
    :param nstereo: max number of generated stereoisomers
    :param nconf: max number of generated conformers
    :param energy: conformers with energy difference from the lowest one greater than the specified value will be discarded.
    :param rms: only conformers with RMS higher then threshold will be kept.
    :param ncpu: number of cpus to use for processing of actives and inactives separately.
    :return:
    """

    comm_path = os.path.join(os.path.dirname(os.path.abspath(in_fname)), 'compounds')
    if not os.path.exists(comm_path):
        os.mkdir(comm_path)

    mol_act = os.path.join(comm_path, 'active.smi')
    mol_inact = os.path.join(comm_path, 'inactive.smi')
    split(in_fname, mol_act, mol_inact)

    procs = []
    for index, fname in enumerate([mol_act, mol_inact]):
        nickname = os.path.splitext(os.path.basename(fname))[0]
        proc = Process(target=common, args=(fname, os.path.join(comm_path, nickname),
                                            nstereo, nconf, energy, rms, rdkit_factory, ncpu, bin_step, nickname))
        procs.append(proc)
        proc.start()

    for proc in procs:
        proc.join()


def entry_point():
    parser = create_parser()
    args = vars(parser.parse_args())
    for o, v in args.items():
        if o == "input": in_fname = v
        if o == "rdkit_factory": rdkit_factory = v if v is not None else None
        if o == "nstereo": nstereo = int(v)
        if o == "nconf": nconf = int(v)
        if o == "energy_cutoff": energy = float(v)
        if o == "rms": rms = float(v) if v is not None else None
        if o == 'bin_step': bin_step = int(v)
        if o == "ncpu": ncpu = int(v)

    main(in_fname=in_fname,
         rdkit_factory=rdkit_factory,
         nstereo=nstereo,
         nconf=nconf,
         energy=energy,
         rms=rms,
         bin_step=bin_step,
         ncpu=ncpu)


if __name__ == '__main__':
    entry_point()
