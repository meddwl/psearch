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

from scripts import gen_stereo_rdkit, gen_conf_rdkit, split
from scripts import create_db


def generate_tautomers(input_fname, output_fname):
    pipe = Popen(
        'cxcalc -g tautomers -f sdf -n true -M true -P true -T true -H 7.4 %s > %s' % (input_fname, output_fname),
        shell=True, stdout=PIPE, stderr=PIPE)
    pipe.wait()


def common(filenames, nconf, energy, rms, gen_tautomers, tolerance, ncpu, fdef_fname, set_name):
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
                               seed=42,
                               verbose=True)

    create_db.main_params(out_fname=None,
                          dbout_fname=filenames[4],
                          smarts_features_fname=None,
                          rdkit_factory=fdef_fname,
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


def main(in_fname, act_threshold, inact_threshold, gen_tautomers,
         nconf, energy, rms, tolerance, fdef_fname, ncpu):

    comm_path = os.path.join(os.path.dirname(os.path.abspath(in_fname)), 'compounds')
    if not os.path.exists(comm_path):
        os.mkdir(comm_path)

    mol_act = os.path.join(comm_path, 'active.smi')
    mol_inact = os.path.join(comm_path, 'inactive.smi')

    split.main(in_fname, mol_act, mol_inact, act_threshold, inact_threshold)

    list_act = [mol_act,
                os.path.join(comm_path, 'active.sdf'),
                os.path.join(comm_path, 'active_stereo.smi'),
                os.path.join(comm_path, 'active_conf.sdf'),
                os.path.join(os.path.split(comm_path)[0], 'active.db')]

    list_inact = [mol_inact,
                  os.path.join(comm_path, 'inactive.sdf'),
                  os.path.join(comm_path, 'inactive_stereo.smi'),
                  os.path.join(comm_path, 'inactive_conf.sdf'),
                  os.path.join(os.path.split(comm_path)[0], 'inactive.db')]

    pa = Process(target=common, args=(list_act, nconf, energy, rms, gen_tautomers, tolerance, ncpu, fdef_fname, 'active'))
    pi = Process(target=common, args=(list_inact, nconf, energy, rms, gen_tautomers, tolerance, ncpu, fdef_fname, 'inactive'))
    pa.start()
    pi.start()
    pa.join()
    pi.join()


if __name__ == '__main__':
    parser = ArgumentParser(description='', formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--in', metavar='input.smi', required=True,
                        help='input file or files.')
    parser.add_argument('-u', '--act_threshold', metavar='act_threshold', default=8.0,
                        help='specify threshold used to determine active compounds.'
                             'Compounds having activity higher or equal to the given'
                             'value will be recognized as active.')
    parser.add_argument('-l', '--inact_threshold', metavar='inact_threshold', default=6.0,
                        help='specify threshold used to determine inactive compounds.'
                             'Compounds having activity less or equal to the given'
                             'value will be recognized as inactive.')
    parser.add_argument('-t', '--gen_tautomers', action='store_true', default=False,
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
    parser.add_argument('-f', '--fdef_fname', metavar='smarts.fdef', default=os.path.join(os.path.dirname(os.path.abspath(__file__)), 'pmapper', 'smarts_features.fdef'),
                        help='fdef-file with pharmacophore feature definition.')
    parser.add_argument('-c', '--ncpu', metavar='cpu_number', default=1,
                        help='number of cpus to use for processing of actives and inactives separately. '
                             'Therefore therefore this number should not exceed max_cpu/2.')

    args = vars(parser.parse_args())
    for o, v in args.items():
        if o == "in": in_fname = v
        if o == "act_threshold": act_threshold = float(v)
        if o == "inact_threshold": inact_threshold = float(v)
        if o == "gen_tautomers": gen_tautomers = v
        if o == "nconf": nconf = int(v)
        if o == "energy_cutoff": energy = float(v)
        if o == "rms": rms = float(v) if v is not None else None
        if o == "tolerance": tolerance = float(v)
        if o == "fdef_fname": fdef_fname = v
        if o == "ncpu": ncpu = int(v)

    main(in_fname=in_fname,
         act_threshold=act_threshold,
         inact_threshold=inact_threshold,
         gen_tautomers=gen_tautomers,
         nconf=nconf,
         energy=energy,
         rms=rms,
         tolerance=tolerance,
         fdef_fname=fdef_fname,
         ncpu=ncpu)