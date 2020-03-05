#!/usr/bin/env python3
# author          : Pavel Polishchuk
# date            : 13.07.16
# license         : BSD-3
#==============================================================================

__author__ = 'Pavel Polishchuk'

import os, time
import sys
import gzip
import shelve
import argparse
from pmapper import utils
from itertools import combinations
from rdkit import Chem
from rdkit.Chem import AllChem
from multiprocessing import Pool, cpu_count
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions
from .read_input import read_input


def prep_input(fname, id_field_name, rdkit_factory, nconf, nstereo, energy, rms, seed, bin_step, box_mol_names):
    input_format = 'smi' if fname is None else None
    for mol, mol_name in read_input(fname, input_format=input_format, id_field_name=id_field_name):
        if mol_name in box_mol_names:
            sys.stderr.write("This {} is meeting the second time in {}. "
                             "This molecule will be omitted.\n".format(mol_name, os.path.basename(fname)))
            continue
        else:
            box_mol_names.append(mol_name)
        yield mol, mol_name, rdkit_factory, nconf, nstereo, energy, rms, seed, bin_step


def map_gen_conf(args):
    return gen_confs(*args)


def remove_confs(mol, energy, rms):

    e = []
    for conf in mol.GetConformers():
        ff = AllChem.MMFFGetMoleculeForceField(mol, AllChem.MMFFGetMoleculeProperties(mol), confId=conf.GetId())
        if ff is None:
            print(Chem.MolToSmiles(mol))
            return
        e.append((conf.GetId(), ff.CalcEnergy()))
    e = sorted(e, key=lambda x: x[1])

    if not e:
        return

    kept_ids = [e[0][0]]
    remove_ids = []

    for item in e[1:]:
        if item[1] - e[0][1] <= energy:
            kept_ids.append(item[0])
        else:
            remove_ids.append(item[0])

    if rms is not None:
        rms_list = [(i1, i2, AllChem.GetConformerRMS(mol, i1, i2)) for i1, i2 in combinations(kept_ids, 2)]
        while any(item[2] < rms for item in rms_list):
            for item in rms_list:
                if item[2] < rms:
                    i = item[1]
                    remove_ids.append(i)
                    break
            rms_list = [item for item in rms_list if item[0] != i and item[1] != i]

    for cid in set(remove_ids):
        mol.RemoveConformer(cid)

    #conformers are reindexed staring with 0 step 1
    for i, conf in enumerate(mol.GetConformers()):
        conf.SetId(i)


def gen_confs(mol, mol_name, rdkit_factory, nconf, nstereo, energy, rms, seed, bin_step):
    set_mol = dict()

    Chem.AssignStereochemistry(mol, flagPossibleStereoCenters=True)
    opts = StereoEnumerationOptions(tryEmbedding=True, maxIsomers=nstereo)
    isomers = tuple(EnumerateStereoisomers(mol, options=opts))

    for i, mol in enumerate(isomers):
        cids = AllChem.EmbedMultipleConfs(mol, numConfs=nconf, maxAttempts=nconf*4, randomSeed=seed)
        for cid in cids:
            AllChem.MMFFOptimizeMolecule(mol, confId=cid)
        remove_confs(mol, energy, rms)

        phs = utils.load_multi_conf_mol(mol, factory=rdkit_factory, bin_step=bin_step)
        set_mol[i] = {
                       'mol': mol,
                       'ph': [ph.get_feature_coords() for ph in phs],
                       'fp': [ph.get_fp() for ph in phs]
                     }
    return mol_name, set_mol


def create_db(in_fname, out_fname, rdkit_factory, id_field_name, nconf, nstereo, energy, rms, ncpu, bin_step, seed, verbose):
    start_time = time.time()
    box_mol_names = []

    output_file_type = None
    if out_fname is not None:

        if os.path.isfile(out_fname):
            os.remove(out_fname)

        if out_fname.lower().endswith('.sdf.gz'):
            writer = gzip.open(out_fname, 'a')
            output_file_type = 'sdf.gz'
        elif out_fname.lower().endswith('.sdf'):
            writer = open(out_fname, 'at')
            output_file_type = 'sdf'
        elif out_fname.lower().endswith('.pkl'):
            writer = open(out_fname, 'wb')
            output_file_type = 'pkl'
        elif out_fname.lower().endswith(''):
            writer = shelve.open(out_fname)
            output_file_type = 'shelve'
            writer['_bin_step'] = bin_step
        else:
            raise Exception("Wrong output file format. Can be only SDF, SDF.GZ or PKL.")

    nprocess = min(cpu_count(), max(ncpu, 1))
    p = Pool(nprocess)

    try:
        for i, (mol_name, set_mol) in enumerate(p.imap_unordered(map_gen_conf, prep_input(in_fname, id_field_name, rdkit_factory,
                                                             nconf, nstereo, energy, rms, seed, bin_step, box_mol_names),
                                                             chunksize=10), 1):
            if output_file_type == 'shelve':
                writer[mol_name] = set_mol

            # elif output_file_type == 'pkl':
            #     pickle.dump((mol, mol_name), writer, -1)
            # else:
            #     mol.SetProp("_Name", mol_name)
            #     string = "$$$$\n".join(Chem.MolToMolBlock(mol, confId=c.GetId()) for c in mol.GetConformers())
            #     if string:   # wrong molecules (no valid conformers) will result in empty string
            #         string += "$$$$\n"
            #         if out_fname is None:
            #             sys.stdout.write(string)
            #             sys.stdout.flush()
            #         else:
            #             writer.write(string.encode("ascii") if output_file_type == 'sdf.gz' else string)
            if verbose and i % 10 == 0:
                sys.stderr.write('\r%i molecules passed/conformers (%is)' % (i, time.time() - start_time))
                sys.stderr.flush()

    finally:
        p.close()

    if out_fname is not None:
        writer.close()
        
    if verbose:
        sys.stderr.write("\n")


def entry_point():
    parser = argparse.ArgumentParser(description='Generate specified number of conformers using RDKit.')
    parser.add_argument('-i', '--in', metavar='input.sdf', required=False, default=None,
                        help='input file with structures to generate conformers. Allowed formats SDF or SMILES. '
                             'if omitted STDIN will be used. STDIN takes only SMILES input (one or two columns).')
    parser.add_argument('-o', '--out', metavar='output.sdf', required=False, default=None,
                        help='output SDF file where conformers are stored. If extension will be SDF.GZ the output file '
                             'will be automatically gzipped. Alternatively for faster storage output can be stored '
                             'in a file with extension PKL. That is pickled storage of tuples (mol, mol_name). '
                             'If the output option will be omitted the output will be done to STDOUT in SDF format.')
    parser.add_argument('-f', '--rdkit_factory', metavar='features.fdef', default=None,
                        help='text file with definition of pharmacophore features in RDKit format. If file name is not '
                             'specified the default file from the script dir will be used. This option has '
                             'a priority over smarts_features.')
    parser.add_argument('-d', '--id_field_name', metavar='field_name', default=None,
                        help='field name of compound ID in input SDF file. If omitted for sdf molecule titles '
                             'will be used or SMILES strings as names.')
    parser.add_argument('-n', '--nconf', metavar='conf_number', default=50,
                        help='number of generated conformers. Default: 50.')
    parser.add_argument('-ns', '--nstereo', metavar='stereo_number', default=5,
                        help='number of generated stereoisomers. Default: 3.')
    parser.add_argument('-e', '--energy_cutoff', metavar='10', default=10,
                        help='conformers with energy difference from the lowest found one higher than the specified '
                             'value will be discarded. Default: 10.')
    parser.add_argument('-r', '--rms', metavar='rms_threshold', default=None,
                        help='only conformers with RMS higher then threshold will be kept. '
                             'Default: None (keep all conformers).')
    parser.add_argument('-b', '--bin_step', default=1,
                        help='binning step. Default: 1.')
    parser.add_argument('-s', '--seed', metavar='random_seed', default=-1,
                        help='integer to init random number generator. Default: -1 (means no seed).')
    parser.add_argument('-c', '--ncpu', metavar='cpu_number', default=1,
                        help='number of cpu to use for calculation. Default: 1.')
    parser.add_argument('-v', '--verbose', action='store_true', default=False,
                        help='print progress to STDERR.')

    args = vars(parser.parse_args())
    for o, v in args.items():
        if o == "in": in_fname = v
        if o == "out": out_fname = v
        if o == "rdkit_factory": rdkit_factory = v
        if o == "id_field_name": id_field_name = v
        if o == "nconf": nconf = int(v)
        if o == "nstereo": nstereo = int(v)
        if o == "ncpu": ncpu = int(v)
        if o == "energy_cutoff": energy = float(v)
        if o == "bin_step": bin_step = float(v)
        if o == "seed": seed = int(v)
        if o == "rms": rms = float(v) if v is not None else None
        if o == "verbose": verbose = v

    create_db(in_fname=in_fname,
              out_fname=out_fname,
              rdkit_factory=rdkit_factory,
              id_field_name=id_field_name,
              nconf=nconf,
              nstereo=nstereo,
              energy=energy,
              rms=rms,
              bin_step=bin_step,
              ncpu=ncpu,
              seed=seed,
              verbose=verbose)


if __name__ == '__main__':
    entry_point()

