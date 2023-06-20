#!/usr/bin/env python3

import os
import sys
import argparse
import pandas as pd
from datetime import datetime
from functools import partial
from itertools import combinations
from multiprocessing import Pool, cpu_count

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions
from pmapper import utils
from psearch.scripts.read_input import read_input
from psearch.database import DB


def create_argparser():
    parser = argparse.ArgumentParser(description='Generates a database of RDKit molecule objects, '
                                                 'coordinates of molecular pharmacophore representations and'
                                                 'pharmacophore fingerprints.',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--input', metavar='FILENAME', required=True, type=str,
                        help='input file of 2D SDF or SMILES format (tab-separated).')
    parser.add_argument('-o', '--db', metavar='FILENAME.dat', required=True, type=str,
                        help='output database file. Should have DAT extension. Database will consist of two files '
                             '.dat and .dir. If there is a database with the same name, then the tool will stop.')
    parser.add_argument('-b', '--bin_step', metavar='NUMERIC', type=int, default=1,
                        help='binning step for pharmacophores creation.')
    parser.add_argument('-s', '--nstereo', metavar='INTEGER', type=int, default=5,
                        help='maximum number of generated stereoisomers per compound (centers with specified '
                             'stereoconfogurations wil not be altered). ')
    parser.add_argument('-n', '--nconf', metavar='INTEGER', type=int, default=50,
                        help='number of generated conformers. ')
    parser.add_argument('-e', '--energy_cutoff', metavar='NUMERIC', type=float, default=None,
                        help='conformers with energy difference from the lowest one greater than the specified '
                             'threshold will be discarded.')
    parser.add_argument('-r', '--rms', metavar='NUMERIC', type=float, default=None,
                        help='only conformers with RMS higher then threshold will be kept. '
                             'Default: None (keep all conformers).')
    parser.add_argument('--seed', metavar='INTEGER', type=int, default=-1,
                        help='integer to init random number generator. Default: -1 (means no seed).')
    parser.add_argument('-p', '--pharm_def', metavar='FILENAME', type=str, default=None,
                        help='pharmacophore feature definition. '
                             'If not specified, default pmapper definitions will be used.')
    parser.add_argument('-c', '--ncpu', metavar='INTEGER', type=int, default=1,
                        help='number of cpu to use for calculation.')
    parser.add_argument('-v', '--verbose', action='store_true', default=False,
                        help='print progress to STDERR.')
    return parser


def check_dupl_input(fname):
    n = 0
    box_mols = []
    for mol, mol_name in read_input(fname):
        if mol in list(zip(*box_mols))[0]:
            n += 1
            sys.stdout.write(
                f'WARNING! The molecule with name {str(mol_name)} meets the second time in the input file. '
                f'This molecule will be omitted\n')
        elif mol_name in list(zip(*box_mols))[1]:
            n += 1
            suffix = 1
            # if there are more than two molecules have the same name
            while mol_name in list(zip(*box_mols))[1]:
                sys.stdout.write(
                    f'WARNING! The molecule ID {str(mol_name)} meets the second time for the distinct molecule SMILES. '
                    f'New molecule ID will be given to this molecule - {str(mol_name)}#{str(suffix)}\n')
                mol_name = f"{str(mol_name)}#{str(suffix)}"
                suffix += 1
            mol.SetProp("_Name", mol_name)
            box_mols.append((mol, mol_name))
        else:
            mol.SetProp("_Name", mol_name)
            box_mols.append((mol, mol_name))
    return list(zip(*box_mols))[0], n


def get_mol(molobj):
    for mol in molobj:
        yield mol


def gen_stereo(mol, num_isomers):
    Chem.AssignStereochemistry(mol, flagPossibleStereoCenters=True)
    opts = StereoEnumerationOptions(tryEmbedding=True, maxIsomers=num_isomers)
    isomers = tuple(EnumerateStereoisomers(mol, options=opts))
    return isomers


def gen_conf(mol, num_confs, seed):
    mol = Chem.AddHs(mol)
    cids = AllChem.EmbedMultipleConfs(mol, numConfs=num_confs, maxAttempts=num_confs*4, randomSeed=seed)
    for cid in cids:
        AllChem.MMFFOptimizeMolecule(mol, confId=cid)
    return mol


def remove_confs(mol, energy, rms):
    e = []
    for conf in mol.GetConformers():
        ff = AllChem.MMFFGetMoleculeForceField(mol, AllChem.MMFFGetMoleculeProperties(mol), confId=conf.GetId())
        if ff is None:
            sys.stderr.write(Chem.MolToSmiles(mol) + ". MMFFGetMoleculeForceField return NONE\n")
            return
        e.append((conf.GetId(), ff.CalcEnergy()))
    e = sorted(e, key=lambda x: x[1])

    if not e:
        return

    kept_ids = [e[0][0]]
    remove_ids = []
    
    if energy is not None:
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
                    remove_ids.append(item[1])
                    rms_list = [i for i in rms_list if i[0] != item[1] and i[1] != item[1]]
                    break

    for cid in set(remove_ids):
        mol.RemoveConformer(cid)

    # conformers are reindexed staring with 0 step 1
    for i, conf in enumerate(mol.GetConformers()):
        conf.SetId(i)


def gen_data(mol, nconf, nstereo, energy, rms, seed, bin_step, pharm_def):
    mol_dict, ph_dict, fp_dict = dict(), dict(), dict()

    isomers = gen_stereo(mol, nstereo)
    for i, mol in enumerate(isomers):
        mol = gen_conf(mol, nconf, seed)
        remove_confs(mol, energy, rms)

        phs = utils.load_multi_conf_mol(mol, smarts_features=pharm_def, bin_step=bin_step)
        mol_dict[i] = mol
        ph_dict[i] = [ph.get_feature_coords() for ph in phs]
        fp_dict[i] = [ph.get_fp() for ph in phs]
    return mol.GetProp("_Name"), mol_dict, ph_dict, fp_dict


def create_db(in_fname, out_fname, nconf, nstereo, energy, rms, ncpu, bin_step, pharm_def, seed, verbose):
    if verbose:
        now = datetime.now()
        date_time = now.strftime("%m/%d/%Y, %H:%M:%S")
        sys.stdout.write(f"Database creation started. {date_time}\n")

    if out_fname.lower().endswith('.dat'):
        db = DB(out_fname, flag='n')
        db.write_bin_step(bin_step)
    else:
        raise Exception("Wrong output file format. Can be only DAT.\n")

    input_molobj, n = check_dupl_input(in_fname)

    nprocess = min(cpu_count(), max(ncpu, 1))
    p = Pool(nprocess)
    try:
        for i, data in enumerate(
                p.imap_unordered(partial(gen_data, nconf, nstereo, energy, rms, seed, bin_step, pharm_def), get_mol(input_molobj),
                                 chunksize=1), 1):
            mol_name, mol_dict, ph_dict, fp_dict = data
            db.write_mol(mol_name, mol_dict)
            db.write_pharm(mol_name, ph_dict)
            db.write_fp(mol_name, fp_dict)

            if i % 200 == 0:
                if verbose:
                    now = datetime.now()
                    date_time = now.strftime("%m/%d/%Y, %H:%M:%S")
                    sys.stdout.write(f"{i} number of molecules were processed. {date_time} \n",)

    finally:
        p.close()

    # create new smi file if the input file has bad molecule structure(-s)
    if n > 0:
        if input_molobj[0].GetProp("Activity"):
            data = [(Chem.MolToSmiles(mol), mol.GetProp("_Name"), mol.GetProp("Activity")) for mol in input_molobj]
            columns = ['smi', 'mol_id', 'activity']
        else:
            data = [(Chem.MolToSmiles(mol), mol.GetProp("_Name")) for mol in input_molobj]
            columns = ['smi', 'mol_id']

        suffix = 1
        new_in_fname = os.path.join(os.path.dirname(in_fname), os.path.basename(in_fname) + f".{str(suffix)}#")
        while os.path.exists(new_in_fname):
            new_in_fname = os.path.join(os.path.dirname(in_fname), os.path.basename(in_fname) + f".{str(suffix)}#")
            suffix += 1
        os.rename(in_fname, new_in_fname)
        df = pd.DataFrame(data=data, columns=columns)
        df.to_csv(os.path.splitext(in_fname)[0] + '.smi', sep='\t', index=None)
        sys.stdout.write(f"\nWARNING! {n} molecules were omitted or renamed comparing with the original input file."
                         f"\nThe original input file will be saved with a new name {os.path.basename(new_in_fname)}."
                         f"\nThe molecules corresponding to the generated database are stored in "
                         f"{os.path.splitext(in_fname)[0]}.smi file\n")


if __name__ == '__main__':
    parser = create_argparser()
    args = parser.parse_args()

    if (args.bin_step < 0) or (args.nstereo <= 0) or (args.nconf <= 0):
        sys.exit("--bin_step, --nstereo, --nconf can not be less 0.\n"
                 "--stereo and/or --nconf can not be set to 0, otherwise, the database will not be created correctly.")

    fdb = os.path.abspath(args.db)
    if os.path.exists(fdb):
        sys.exit(f"Database with this {fdb} name already exists")
    else:
        os.makedirs(os.path.dirname(fdb), exist_ok=True)

    create_db(in_fname=os.path.abspath(args.input),
              out_fname=fdb,
              nconf=args.nconf,
              nstereo=args.nstereo,
              energy=args.energy_cutoff,
              rms=args.rms,
              bin_step=args.bin_step,
              pharm_def=args.pharm_def,
              ncpu=args.ncpu,
              seed=args.seed,
              verbose=args.verbose)
