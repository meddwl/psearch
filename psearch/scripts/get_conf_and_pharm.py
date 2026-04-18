#!/usr/bin/env python3
# author          : Alina Kutlushina
# date            : 10.06.2022
# license         : BSD-3
# ==============================================================================

__author__ = 'Alina Kutlushina'

import os
from rdkit import Chem
# from psearch.database import DB
from psearch.database import load_model
import argparse


class RawTextArgumentDefaultsHelpFormatter(argparse.RawTextHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
    """Combined argparse formatter that preserves raw text line breaks and appends default values to help strings."""
    pass


def extractor(db_path, mol_ids, stereo_ids, conf_ids):
    """Retrieve specific conformers and their pharmacophore features from the database.

    Args:
        db_path: Path to the psearch SQLite database file.
        mol_ids: List of compound identifier strings.
        stereo_ids: List of stereo_id integers (one per entry in mol_ids).
        conf_ids: List of conformer_id integers (one per entry in mol_ids).

    Yields:
        Tuples of (rdkit.Mol, conformer_name, conf_id, pharm_string) where
        conformer_name is '{mol_id}-s{stereo_id}-c{conf_id}' and pharm_string
        is the pharmacophore features formatted for .xyz output.
    """
    for mol, stereo, conf in zip(mol_ids, stereo_ids, conf_ids):
        # record = db.get_mol(mol), db.get_pharm(mol) replaced by load_model
        record = load_model(db_path, mol)
        cmp = record.mol_dict[stereo]
        mname = f"{mol}-s{stereo}-c{conf}"
        cmp.SetProp("_Name", mname)

        # pharm = db.get_pharm(mol)[stereo][conf]
        pharm = record.pharm_dict[stereo][conf]
        pharm = "\n\n" + "\n".join([i[0] + ' ' + ' '.join(map(str, i[1])) for i in pharm])
        yield cmp, mname, conf, pharm


def main():
    parser = argparse.ArgumentParser(description="""
    Extract a particular conformer of a molecule in sdf format file and its pharmacophore in .xyz format file 
    from a psearch database for a required conformer(-s) by giving molecule id(-s), stereo id(-s) and conformer id(-s). 
             """, formatter_class=RawTextArgumentDefaultsHelpFormatter)
    parser.add_argument('-d', '--dbdir', metavar='database.db', required=True, type=str,
                        help='path to a psearch database')
    parser.add_argument('-m', '--mol_id', metavar='molecule_name', nargs='+', required=True, type=str,
                        help='molecule ID of a required conformer(-s) in the psearch database')
    parser.add_argument('-s', '--stereo_id', metavar='stereo_id', nargs='+', required=True, type=int,
                        help='stereo ID of a required conformer(-s) in the psearch database')
    parser.add_argument('-c', '--conf_id', metavar='conformer_id', nargs='+', required=True, type=int,
                        help='conformer ID of a required conformer(-s) in the psearch database')
    parser.add_argument('-o', '--output', metavar='DIRNAME', required=True, type=str,
                        help='a folder path where will be saved the results')

    args = parser.parse_args()

    # db = DB(os.path.abspath(args.dbdir))
    db_path = os.path.abspath(args.dbdir)
    for cmp, cmp_name, conf_id, pharm in extractor(db_path, args.mol_id, args.stereo_id, args.conf_id):
        writer = Chem.SDWriter(os.path.join(os.path.abspath(args.output), cmp_name + '.sdf'))
        writer.write(cmp, confId=conf_id)

        with open(os.path.join(os.path.abspath(args.output), cmp_name), 'a') as f:
            f.write(pharm)


if __name__ == '__main__':
    main()
