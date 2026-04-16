#!/usr/bin/env python3
# author          : Pavel Polishchuk
# date            : 04.03.20
# license         : BSD-3
# ==============================================================================

__author__ = 'Pavel Polishchuk'

import os
import shelve


class DB:
    """Context-manager wrapper around a shelve database storing molecular conformers,
    pharmacophore feature sets, and binary fingerprints for a compound library.
    """

    def __init__(self, fname, flag='c'):
        """Open (or create) the psearch shelve database.

        Args:
            fname: Path to the database file (the extension is stripped automatically;
                   shelve creates .dat and .dir sidecar files).
            flag: shelve open flag — 'c' creates or opens for read/write (default),
                  'r' opens read-only, 'n' always creates a new empty database.
        """
        self.__db = shelve.open(os.path.splitext(fname)[0], flag=flag, protocol=4)

    def __enter__(self):
        """Return self to support use as a context manager."""
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        """Close the underlying shelve database on context exit."""
        self.__db.close()

    def write_bin_step(self, bin_step):
        """Persist the binning step (Å) used to discretise pharmacophore coordinates."""
        self.__db['_bin_step'] = bin_step

    def write_mol(self, mol_name, mol_dict):
        """Store 3D conformer objects for a compound.

        Args:
            mol_name: Compound identifier string.
            mol_dict: Dict mapping stereo_id (int) to an RDKit Mol with embedded conformers.
        """
        self.__db[f'{mol_name}_mol'] = mol_dict

    def write_pharm(self, mol_name, pharm_dict):
        """Store pharmacophore feature coordinates for a compound.

        Args:
            mol_name: Compound identifier string.
            pharm_dict: Dict mapping stereo_id (int) to a list of feature-coordinate tuples,
                        one list per conformer: [[(label, (x, y, z)), …], …].
        """
        self.__db[f'{mol_name}_pharm'] = pharm_dict

    def write_fp(self, mol_name, fp_dict):
        """Store binary pharmacophore fingerprints for a compound.

        Args:
            mol_name: Compound identifier string.
            fp_dict: Dict mapping stereo_id (int) to a list of fingerprint dicts,
                     one dict per conformer.
        """
        self.__db[f'{mol_name}_fp'] = fp_dict

    def get_bin_step(self):
        """Return the binning step (Å) stored in the database."""
        return self.__db['_bin_step']

    def get_mol(self, mol_name):
        """Return the stored RDKit Mol objects for mol_name, keyed by stereo_id (int)."""
        return self.__db[f'{mol_name}_mol']

    def get_pharm(self, mol_name):
        """Return pharmacophore feature coordinate data for mol_name, keyed by stereo_id (int)."""
        return self.__db[f'{mol_name}_pharm']

    def get_fp(self, mol_name):
        """Return binary pharmacophore fingerprint data for mol_name, keyed by stereo_id (int)."""
        return self.__db[f'{mol_name}_fp']

    def get_mol_names(self):
        """Return a tuple of all compound identifiers stored in the database."""
        names = list(self.__db.keys())
        names.remove('_bin_step')
        names = [n[:-4] for n in names if n.endswith('_mol')]
        return tuple(names)

    def get_conf_count(self, mol_name):
        """Return the total number of conformers (across all stereoisomers) stored for mol_name."""
        fps = self.get_fp(mol_name)
        return sum(len(item) for item in fps.values())
