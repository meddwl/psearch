#!/usr/bin/env python3
# author          : Alina Kutlushina
# date            : 01.05.2018
# license         : BSD-3
#==============================================================================

import os
import sys
import time

from subprocess import Popen, PIPE
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from multiprocessing import Process

import gen_stereo_rdkit, gen_conf_rdkit, split
import create_db


def generate_tautomers(input_fname, output_fname):
    pipe = Popen(
        'cxcalc -g tautomers -f sdf -n true -M true -P true -T true -H 7.4 %s > %s' % (input_fname, output_fname),
        shell=True, stdout=PIPE, stderr=PIPE)
    pipe.wait()


def common(filenames, nconf, energy, rms, rdkit_factory, gen_tautomers, tolerance, ncpu, set_name):
    start = time.time()

    if gen_tautomers:
        generate_tautomers(filenames[0], filenames[1])
        input_fname = filenames[1]
    else:
        input_fname = filenames[0]

    gen_stereo_rdkit.main_params(in_fname=input_fname,
                                 out_fname=filenames[2],
                                 tetrahedral=True,
                                 double_bond=True,
                                 max_undef=-1,
                                 id_field_name=None,
                                 ncpu=ncpu,
                                 verbose=True)

    gen_conf_rdkit.main_params(in_fname=filenames[2],
                               out_fname=filenames[3],
                               id_field_name=None,
                               nconf=nconf,
                               energy=energy,
                               rms=rms,
                               ncpu=ncpu,
                               seed=-1,
                               verbose=True)

    create_db.main_params(out_fname=None,
                          dbout_fname=filenames[4],
                          smarts_features_fname=None, 
                          rdkit_factory=rdkit_factory,
                          conformers_fname=filenames[3],
                          bin_step=1,
                          rewrite_db=True,
                          store_coords=True,
                          id_field_name=None,
                          stereo_id=True,
                          fp=True,
                          nohash=True,
                          verbose=True,
                          ncpu=ncpu,
                          tolerance=tolerance)

    sys.stderr.write('prepare {} dataset ({}s)'.format(set_name, time.time() - start))


def main(in_fname, split_dataset, rdkit_factory, gen_tautomers, nconf, energy, rms, tolerance, ncpu):
    """
    launches the entire cycle of data preprocessing: generation of stereoisomers, conformers and a database
    :param in_fname: input .smi file containing information about SMILES, compounds id and its activity status
    :param split_dataset: if True will splited input dasets into active and inactive sets else will not
    :param rdkit_factory: text file with definition of pharmacophore features in RDKit format.
    :param gen_tautomers: if True will generate tautomers for every compounds else will not
    :param nconf: max number of generated conformers
    :param energy: conformers with energy difference from the lowest one greater than the specified value will be discarded.
    :param rms: only conformers with RMS higher then threshold will be kept.
    :param tolerance: tolerance volume for the calculation of the stereo sign. (for more information read helper)
    :param ncpu: number of cpus to use for processing of actives and inactives separately.
    :return:
    """

    comm_path = os.path.join(os.path.dirname(os.path.abspath(in_fname[0])), 'compounds')
    if not os.path.exists(comm_path):
        os.mkdir(comm_path)

    if split_dataset:
        mol_act = os.path.join(comm_path, 'active.smi')
        mol_inact = os.path.join(comm_path, 'inactive.smi')
        split.main(in_fname[0], mol_act, mol_inact)
        #in_fname = [mol_act, mol_inact]
        in_fname = [mol_act]

    procs = []
    for index, fname in enumerate(in_fname):
        nickname = os.path.basename(fname).split('.')[0]
        list_ts = [fname,
                   os.path.join(comm_path, '{}_taut.sdf'.format(nickname)),
                   os.path.join(comm_path, '{}_stereo.smi'.format(nickname)),
                   os.path.join(comm_path, '{}_conf.sdf'.format(nickname)),
                   os.path.join(comm_path, '{}.db'.format(nickname))]
        proc = Process(target=common, args=(list_ts, nconf, energy, rms, rdkit_factory,
                                            gen_tautomers, tolerance, ncpu, nickname))
        procs.append(proc)
        proc.start()

    for proc in procs:
        proc.join()


if __name__ == '__main__':
    parser = ArgumentParser(description='', formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--input', metavar='input.smi', nargs='+', type=str, required=True,
                        help='input smi file or multiple files')
    parser.add_argument('-s', '--split_input_file', action='store_true', default=True,
                        help='if True will splited input dasets into active and inactive sets if False will not')
    parser.add_argument('-f', '--rdkit_factory', metavar='features.fdef', default=None,
                        help='text file with definition of pharmacophore features in RDKit format. If file name is not '
                             'specified the default file from the script dir will be used. This option has '
                             'a priority over smarts_features.')
    parser.add_argument('-g', '--gen_tautomers', action='store_true', default=False,
                        help='if True tautomers at pH 7.4 are generated using Chemaxon.')
    parser.add_argument('-n', '--nconf', metavar='conf_number', default=100,
                        help='number of generated conformers.')
    parser.add_argument('-e', '--energy_cutoff', metavar='100', default=100,
                        help='conformers with energy difference from the lowest one greater than the specified '
                             'value will be discarded.')
    parser.add_argument('-r', '--rms', metavar='rms_threshold', default=0.5,
                        help='only conformers with RMS higher then threshold will be kept.')
    parser.add_argument('-tol', '--tolerance', default=0,
                        help='tolerance volume for the calculation of the stereo sign. If the volume of the '
                             'tetrahedron created by four points less than tolerance then those points are considered '
                             'lying on the same plane (flat; stereo sign is 0).')
    parser.add_argument('-c', '--ncpu', metavar='cpu_number', default=1,
                        help='number of cpus to use for processing of actives and inactives separately. ')

    args = vars(parser.parse_args())
    for o, v in args.items():
        if o == "input": in_fname = v
        if o == "split_input_file": split_dataset = v
        if o == "rdkit_factory": rdkit_factory = v
        if o == "gen_tautomers": gen_tautomers = v
        if o == "nconf": nconf = int(v)
        if o == "energy_cutoff": energy = float(v)
        if o == "rms": rms = float(v) if v is not None else None
        if o == "tolerance": tolerance = float(v)
        if o == "ncpu": ncpu = int(v)


    main(in_fname=in_fname,
         split_dataset=split_dataset,
         rdkit_factory=rdkit_factory,
         gen_tautomers=gen_tautomers,
         nconf=nconf,
         energy=energy,
         rms=rms,
         tolerance=tolerance,
         ncpu=ncpu)
