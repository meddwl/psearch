#!/usr/bin/env python3
# author          : Pavel Polishchuk
# date            : 04.03.20
# license         : BSD-3
# ==============================================================================

__author__ = 'Pavel Polishchuk'

import os
import shelve


class DB:
    def __init__(self, fname):
        self.__db = shelve.open(os.path.splitext(fname)[0])

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.__db.close()

    def write_bin_step(self, bin_step):
        self.__db['_bin_step'] = bin_step

    def write_mol(self, mol_name, mol_dict):
        self.__db[f'{mol_name}_mol'] = mol_dict

    def write_pharm(self, mol_name, pharm_dict):
        self.__db[f'{mol_name}_pharm'] = pharm_dict

    def write_fp(self, mol_name, fp_dict):
        self.__db[f'{mol_name}_fp'] = fp_dict

    def get_bin_step(self):
        return self.__db['_bin_step']

    def get_mol(self, mol_name):
        return self.__db[f'{mol_name}_mol']

    def get_pharm(self, mol_name):
        return self.__db[f'{mol_name}_pharm']

    def get_fp(self, mol_name):
        return self.__db[f'{mol_name}_fp']

    def get_mol_names(self):
        names = list(self.__db.keys())
        names.remove('_bin_step')
        names = [n[:-4] for n in names if n.endswith('_mol')]
        return tuple(names)
