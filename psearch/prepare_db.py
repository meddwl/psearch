#!/usr/bin/env python3
# author          : Alina Kutlushina
# date            : 01.05.2018
# license         : BSD-3
#==============================================================================
import os
import sys

from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

from .scripts import gen_stereo_rdkit, gen_conf_rdkit, split
from .scripts import create_db


def create_parser():
    parser = ArgumentParser(description='Create database for virtual screening from input structures. '
                                        'Stereoisomers and conformers will be generated during processing. '
                                        'Intermediate files will be stored in the same directory as output database.',
                            formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--input', metavar='input.smi', type=os.path.abspath, required=True,
                        help='input file of SMI or SDF formats, recognized by extension.')
    parser.add_argument('-o', '--output', metavar='output.db', type=os.path.abspath, required=True,
                        help='input file of SMI or SDF formats, recognized by extension.')
    parser.add_argument('-u', '--max_undef', metavar='INTEGER', required=False, default=-1, type=int,
                        help='maximum allowed number of unspecified stereocenters and/or double bonds. '
                             'if compound contains greater number of them it will be discarded. '
                             'Default: -1, all possible stereoisomers will be enumerated '
                             '(beware of combinatorial explosion).')
    parser.add_argument('-n', '--nconf', metavar='INTEGER', required=False, default=100, type=int,
                        help='number of generated conformers.')
    parser.add_argument('-e', '--energy_cutoff', metavar='NUMERIC', required=False, default=100, type=float,
                        help='conformers with energy difference from the lowest one greater than the specified '
                             'value will be discarded.')
    parser.add_argument('-r', '--rms', metavar='NUMERIC', required=False, default=0.5, type=float,
                        help='only conformers with RMS higher then threshold will be kept.')
    parser.add_argument('-s', '--smarts_features', metavar='smarts.txt', required=False, default=None, type=str,
                        help='text file with definition of pharmacophore features. Tab-separated two columns: '
                             'first - SMARTS, second - feature name. If file name is not specified the default SMARTS '
                             'patterns will be used.')
    parser.add_argument('-f', '--rdkit_factory', metavar='features.fdef', required=False, default=None, type=str,
                        help='text file with definition of pharmacophore features in RDKit format. If file name is not '
                             'specified the default SMARTS patterns will be used. This option has a priority over '
                             'smarts_features if both were specified.')
    parser.add_argument('-c', '--ncpu', metavar='INTEGER', required=False, default=1, type=int,
                        help='number of cpus to use for processing of actives and inactives separately. ')
    parser.add_argument('-v', '--verbose', required=False, action='store_true', default=False,
                        help='print progress.')
    return parser


def entry_point():
    parser = create_parser()
    args = parser.parse_args()

    if not os.path.isdir(os.path.dirname(args.output)):
        os.makedirs(os.path.dirname(args.output), exist_ok=True)

    if args.verbose:
        sys.stderr.write('Stereoisomers enumeration started\n')
        sys.stderr.flush()

    if args.max_undef != 0:
        fname_stereo = os.path.splitext(args.output)[0] + '_stereo.smi'
        gen_stereo_rdkit.main_params(in_fname=args.input,
                                     out_fname=fname_stereo,
                                     tetrahedral=True,
                                     double_bond=True,
                                     max_undef=args.max_undef,
                                     id_field_name=None,
                                     ncpu=args.ncpu,
                                     verbose=args.verbose)
    else:
        fname_stereo = args.input

    if args.verbose:
        sys.stderr.write('Conformers generation started\n')
        sys.stderr.flush()

    fname_conf = os.path.splitext(args.output)[0] + '_conf.sdf'

    gen_conf_rdkit.main_params(in_fname=fname_stereo,
                               out_fname=fname_conf,
                               id_field_name=None,
                               nconf=args.nconf,
                               energy=args.energy_cutoff,
                               rms=args.rms,
                               ncpu=args.ncpu,
                               seed=-1,
                               verbose=args.verbose)

    if args.verbose:
        sys.stderr.write('Database creation started\n')
        sys.stderr.flush()

    create_db.main_params(dbout_fname=args.output,
                          smarts_features_fname=args.smarts_features,
                          rdkit_factory=args.rdkit_factory,
                          conformers_fname=fname_conf,
                          bin_step=1,
                          rewrite_db=True,
                          id_field_name=None,
                          stereo_id=True,
                          verbose=args.verbose,
                          ncpu=args.ncpu)


if __name__ == '__main__':
    entry_point()
