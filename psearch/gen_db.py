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
    """Build the CLI argument parser for gen_db."""
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
                        help='Bin width (Å) for discretising pharmacophore feature coordinates. '
                             'Smaller values increase fingerprint resolution. Default: 1 Å.')
    parser.add_argument('-s', '--nstereo', metavar='INTEGER', type=int, default=5,
                        help='Maximum number of stereoisomers to generate per compound. '
                             'Stereocentres with explicitly specified configurations will not be altered.')
    parser.add_argument('-n', '--nconf', metavar='INTEGER', type=int, default=50,
                        help='number of generated conformers. ')
    parser.add_argument('-e', '--energy_cutoff', metavar='NUMERIC', type=float, default=None,
                        help='conformers with energy difference from the lowest one greater than the specified '
                             'threshold will be discarded.')
    parser.add_argument('-r', '--rms', metavar='NUMERIC', type=float, default=None,
                        help='Only conformers with pairwise RMS higher than this threshold (Å) will be kept. '
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
    """Check the input file for duplicate SMILES or molecule IDs.

    Duplicate SMILES are skipped and flagged; duplicate IDs for distinct SMILES are
    renamed by appending a '#N' suffix. Returns parallel lists used to build the database.

    Args:
        fname: Path to the 2D SDF or tab-separated SMILES input file.

    Returns:
        Tuple of (mols, smis, mol_ids, flags) where flags[i] is True if molecule i
        was a duplicate or was renamed.
    """
    n = 0
    box_mols, box_smi, box_mol_ids, box_flags = [], [], [], []
    smi_set = set()   # O(1) lookup instead of O(n) list search
    mol_id_set = set()
    for mol, mol_name in read_input(fname):
        canon_smi = Chem.MolToSmiles(mol)
        if canon_smi in smi_set:
            n += 1
            sys.stdout.write(
                f'WARNING! The molecule with name {str(mol_name)} meets the second time in the input file. '
                f'This molecule will be omitted\n')
            box_mols.append("duplicate")
            box_smi.append("duplicate")
            box_mol_ids.append(mol_name)
            box_flags.append(True)
        elif mol_name in mol_id_set:
            n += 1
            suffix = 1
            # if there are more than two molecules have the same name
            while mol_name in mol_id_set:
                sys.stdout.write(
                    f'WARNING! The molecule ID {str(mol_name)} meets the second time for the distinct molecule SMILES. '
                    f'New molecule ID will be given to this molecule - {str(mol_name.split("#")[0])}#{str(suffix)}\n')
                mol_name = f"{str(mol_name.split('#')[0])}#{str(suffix)}"
                suffix += 1
            box_mols.append(mol)
            box_smi.append(canon_smi)
            box_mol_ids.append(mol_name)
            box_flags.append(True)
            smi_set.add(canon_smi)
            mol_id_set.add(mol_name)
        else:
            box_mols.append(mol)
            box_smi.append(canon_smi)
            box_mol_ids.append(mol_name)
            box_flags.append(False)
            smi_set.add(canon_smi)
            mol_id_set.add(mol_name)

    return box_mols, box_smi, box_mol_ids, box_flags


def get_mol(mols):
    """Yield (mol, mol_name) tuples from mols, skipping entries marked as duplicates."""
    for mol, mol_name in mols:
        if mol == 'duplicate':
            continue
        yield mol, mol_name


def gen_stereo(mol, num_isomers):
    """Enumerate stereoisomers of mol, leaving explicitly defined stereocentres unchanged.

    Args:
        mol: RDKit Mol object (2D, no explicit stereo required).
        num_isomers: Maximum number of stereoisomers to generate.

    Returns:
        Tuple of RDKit Mol objects, one per enumerated stereoisomer.
    """
    Chem.AssignStereochemistry(mol, flagPossibleStereoCenters=True)
    opts = StereoEnumerationOptions(tryEmbedding=True, maxIsomers=num_isomers)
    isomers = tuple(EnumerateStereoisomers(mol, options=opts))
    return isomers


def gen_conf(mol, num_confs, seed):
    """Generate 3D conformers for mol using ETKDGv3 followed by MMFF optimisation.

    Args:
        mol: RDKit Mol object (typically a stereoisomer from gen_stereo).
        num_confs: Number of conformers to embed.
        seed: Random seed for ETKDGv3 (use -1 for no fixed seed).

    Returns:
        RDKit Mol with hydrogens added and conformers embedded and minimised.
    """
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    params.maxAttempts = num_confs * 4
    AllChem.EmbedMultipleConfs(mol, numConfs=num_confs, params=params)
    AllChem.MMFFOptimizeMoleculeConfs(mol, numThreads=0)  # numThreads=0 uses all available cores
    return mol


def remove_confs(mol, energy, rms):
    """Remove conformers that fail energy or RMS diversity filters, in-place.

    Conformers with MMFF energy more than `energy` kcal/mol above the lowest-energy
    conformer are removed first. Then, among the surviving conformers, redundant ones
    with pairwise RMS below `rms` Å are removed (keeping the lower-energy member of
    each close pair). Conformers are reindexed starting from 0 after pruning.

    Args:
        mol: RDKit Mol with embedded and MMFF-minimised conformers.
        energy: Maximum allowed energy difference (kcal/mol) from the lowest conformer.
                None disables the energy filter.
        rms: Minimum pairwise RMS (Å) required to retain two conformers.
             None disables the RMS filter.
    """
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


def gen_data(mol_cid, nconf, nstereo, energy, rms, seed, bin_step, pharm_def):
    """Generate stereoisomers, embed conformers, and compute pharmacophores for one compound.

    Enumerates up to `nstereo` stereoisomers, embeds up to `nconf` conformers per
    stereoisomer (filtered by `energy` and `rms`), then computes pmapper pharmacophore
    feature coordinates and binary fingerprints for each conformer.

    Args:
        mol_cid: Tuple of (rdkit.Mol, mol_name) for the input 2D compound.
        nconf: Number of conformers to generate per stereoisomer.
        nstereo: Maximum number of stereoisomers to enumerate.
        energy: MMFF energy filter cutoff (kcal/mol); None to disable.
        rms: RMS diversity filter (Å); None to disable.
        seed: Random seed for conformer embedding (-1 for no seed).
        bin_step: Bin width (Å) for pharmacophore coordinate discretisation.
        pharm_def: Path to a custom pmapper feature-definition file, or None to use defaults.

    Returns:
        Tuple of (mol_name, mol_dict, ph_dict, fp_dict) keyed by stereo_id.
    """
    mol_dict, ph_dict, fp_dict = dict(), dict(), dict()
    mol, mol_name = mol_cid

    isomers = gen_stereo(mol, nstereo)
    for i, mol in enumerate(isomers):
        mol = gen_conf(mol, nconf, seed)
        remove_confs(mol, energy, rms) # what if it returns None?

        phs = utils.load_multi_conf_mol(mol, smarts_features=pharm_def, bin_step=bin_step)
        mol_dict[i] = mol
        ph_dict[i] = [ph.get_feature_coords() for ph in phs]
        fp_dict[i] = [ph.get_fp() for ph in phs]
    return mol_name, mol_dict, ph_dict, fp_dict


def create_db(in_fname, out_fname, nconf, nstereo, energy, rms, ncpu, bin_step, pharm_def, seed, verbose):
    """Orchestrate parallel generation of conformers and pharmacophore fingerprints.

    Reads all compounds from `in_fname`, checks for duplicate SMILES/IDs, then
    processes each compound in parallel (gen_data) and writes results to the psearch
    shelve database at `out_fname`. If any molecules were renamed or skipped, a
    corrected SMILES file is written alongside the input file.

    Args:
        in_fname: Path to the 2D SDF or tab-separated SMILES input file.
        out_fname: Path to the output database file (must have .dat extension).
        nconf: Number of conformers to generate per stereoisomer.
        nstereo: Maximum number of stereoisomers per compound.
        energy: MMFF energy filter cutoff (kcal/mol); None to disable.
        rms: RMS diversity filter (Å); None to disable.
        ncpu: Number of parallel worker processes.
        bin_step: Bin width (Å) for pharmacophore discretisation.
        pharm_def: Path to a custom pmapper feature-definition file, or None for defaults.
        seed: Random seed for conformer embedding (-1 for no seed).
        verbose: If True, print progress and timing information to stdout.
    """
    if verbose:
        now = datetime.now()
        date_time = now.strftime("%m/%d/%Y, %H:%M:%S")
        sys.stdout.write(f"Database creation started. {date_time}\n")

    if out_fname.lower().endswith('.dat'):
        db = DB(out_fname, flag='n')
        db.write_bin_step(bin_step)
    else:
        raise Exception("Wrong output file format. Can be only DAT.\n")

    mols, smis, cids, flags = check_dupl_input(in_fname)
    nprocess = min(cpu_count(), max(ncpu, 1))
    try:
        with Pool(nprocess) as p:
            for i, data in enumerate(
                    p.imap_unordered(partial(gen_data, nconf=nconf, nstereo=nstereo, energy=energy, rms=rms, seed=seed,
                                             bin_step=bin_step, pharm_def=pharm_def), get_mol(zip(mols, cids)), chunksize=1), 1):
                if not data:
                    continue
                mol_name, mol_dict, ph_dict, fp_dict = data
                db.write_mol(mol_name, mol_dict)
                db.write_pharm(mol_name, ph_dict)
                db.write_fp(mol_name, fp_dict)

                if i % 200 == 0:
                    if verbose:
                        now = datetime.now()
                        date_time = now.strftime("%m/%d/%Y, %H:%M:%S")
                        sys.stdout.write(f"{i} number of molecules were processed. {date_time} \n")

        # create new smi file if the input file has bad molecule structure(-s)
        if sum(flags) > 0:
            sys.stdout.write(
                f"\nWARNING! {str(sum(flags))} molecules were omitted and/or renamed "
                f"comparing with the original input file.\n")

            suffix = 2
            new_in_fname = f"{os.path.splitext(in_fname)[0]}-updated.smi"
            while os.path.exists(new_in_fname):
                new_in_fname = f"{os.path.splitext(in_fname)[0]}-updated{str(suffix)}.smi"
                suffix += 1
            df_cid = pd.DataFrame(data={'smi': smis, 'cid': cids, 'if_changed': flags})
            df_cid.to_csv(new_in_fname, sep='\t', index=None)

            sys.stdout.write(
                f"The molecules corresponding to the generated database are stored in {new_in_fname} file\n\n")
    finally:
        if verbose:
            now = datetime.now()
            date_time = now.strftime("%m/%d/%Y, %H:%M:%S")
            sys.stdout.write(f"Database is created. {date_time}\n")


def entry_point():
    """CLI entry point for gen_db: parse arguments and call create_db."""
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


if __name__ == '__main__':
    entry_point()
