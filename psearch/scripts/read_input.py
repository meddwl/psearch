#!/usr/bin/env python3
# author          : Pavel Polishchuk
# date            : 13.07.16
# license         : BSD-3
#==============================================================================

__author__ = 'pavel'

import os
import sys
import gzip
import pickle
from rdkit import Chem
from rdkit.Chem.PropertyMol import PropertyMol
from io import BytesIO


def __read_pkl(fname):
    """Yield (rdkit.Mol, mol_name) tuples from a pickled molecule file."""
    with open(fname, 'rb') as f:
        while True:
            try:
                yield pickle.load(f)
            except EOFError:
                break


def __read_sdf(fname, input_format, id_field_name=None, sanitize=True, removeHs=True):
    """Yield (rdkit.Mol, mol_name) tuples from an SDF or gzipped SDF file.

    The molecule name is taken from the SD property specified by id_field_name, falling
    back to the molecule title (_Name), or to the canonical SMILES if the title is empty.
    """
    if input_format == 'sdf':
        suppl = Chem.SDMolSupplier(fname, sanitize=sanitize, removeHs=removeHs)
    elif input_format == 'sdf.gz':
        suppl = Chem.ForwardSDMolSupplier(gzip.open(fname), sanitize=sanitize, removeHs=removeHs)
    else:
        return
    for mol in suppl:
        if mol is not None:
            if id_field_name is not None:
                mol_title = mol.GetProp(id_field_name)
            else:
                if mol.GetProp("_Name"):
                    mol_title = mol.GetProp("_Name")
                else:
                    mol_title = Chem.MolToSmiles(mol, isomericSmiles=True)
            yield PropertyMol(mol), mol_title


def __read_smiles(fname, sanitize=True):
    """Yield (rdkit.Mol, mol_name) tuples from a tab-separated SMILES file (columns: SMILES, mol_id)."""
    with open(fname) as f:
        for line in f:
            tmp = line.strip().split()
            mol = Chem.MolFromSmiles(tmp[0], sanitize=sanitize)
            if mol is not None:
                if len(tmp) > 1:
                    mol_title = tmp[1]
                else:
                    mol_title = Chem.MolToSmiles(mol, isomericSmiles=True)
                yield mol, mol_title


def __read_stdin_smiles(sanitize=True):
    """Yield (rdkit.Mol, mol_name) tuples by reading tab-separated SMILES from standard input."""
    line = sys.stdin.readline()
    while line:
        tmp = line.strip().split()
        mol = Chem.MolFromSmiles(tmp[0], sanitize=sanitize)
        if mol is not None:
                if len(tmp) > 1:
                    mol_title = tmp[1]
                else:
                    mol_title = Chem.MolToSmiles(mol, isomericSmiles=True)
                yield mol, mol_title
        line = sys.stdin.readline()


def __read_stdin_sdf(sanitize=True, removeHs=True):
    """Yield (rdkit.Mol, mol_name) tuples by reading an SDF stream from standard input."""
    molblock = ''
    line = sys.stdin.readline()
    while line:
        molblock += line
        if line == '$$$$\n':
            mol = [x for x in Chem.ForwardSDMolSupplier(BytesIO(molblock.encode('utf-8')), sanitize=sanitize, removeHs=removeHs)][0]
            mol_title = molblock.split('\n', 1)[0]
            if not mol_title:
                mol_title = Chem.MolToSmiles(mol, isomericSmiles=True)
            yield mol, mol_title
            molblock = ''
        line = sys.stdin.readline()


# def read_input(fname, id_field_name=None, stdin_format=None, sanitize=True):
#     if fname is None:
#         if stdin_format == 'smi':
#             suppl = read_stdin_smiles()
#         elif stdin_format == 'sdf':
#             suppl = read_stdin_sdf(sanitize=sanitize)
#         else:
#             raise Exception("Cannot read STDIN. STDIN format should be specified explicitly: smi or sdf.")
#     elif fname.lower().endswith('.sdf') or fname.lower().endswith('.sdf.gz'):
#         suppl = read_sdf(os.path.abspath(fname), id_field_name=id_field_name, sanitize=sanitize)
#     elif fname.lower().endswith('.smi') or fname.lower().endswith('.smiles'):
#         suppl = read_smiles(os.path.abspath(fname))
#     elif fname.lower().endswith('.pkl'):
#         suppl = read_pkl(os.path.abspath(fname))
#     else:
#         raise Exception("File extension can be only SDF, SMI or SMILES")
#     for mol, mol_name in suppl:
#         yield mol, mol_name


def read_input(fname, input_format=None, id_field_name=None, sanitize=True, removeHs=True):
    """Yield (rdkit.Mol, mol_name) tuples from a molecule file or standard input.

    Args:
        fname: Path to the input file, or None to read from STDIN.
        input_format: File format string ('sdf', 'sdf.gz', 'smi', 'smiles', or 'pkl').
                      Inferred from the file extension when fname is not None.
                      Must be provided explicitly when reading from STDIN.
        id_field_name: Name of the SD property to use as the molecule identifier (SDF only).
                       Falls back to the molecule title (_Name) if None.
        sanitize: If True (default), RDKit sanitizes each molecule on read.
        removeHs: If True (default), explicit hydrogens are removed on read (SDF only).

    Yields:
        Tuples of (rdkit.Mol, mol_name str).
    """
    if input_format is None:
        tmp = os.path.basename(fname).split('.')
        if tmp[-1] == 'gz':
            input_format = '.'.join(tmp[-2:])
        else:
            input_format = tmp[-1]
    input_format = input_format.lower()
    if fname is None:    # handle STDIN
        if input_format == 'sdf':
            suppl = __read_stdin_sdf(sanitize=sanitize, removeHs=removeHs)
        elif input_format == 'smi':
            suppl = __read_stdin_smiles(sanitize=sanitize)
        else:
            raise Exception("Input STDIN format '%s' is not supported. It can be only sdf, smi." % input_format)
    elif input_format in ("sdf", "sdf.gz"):
        suppl = __read_sdf(os.path.abspath(fname), input_format, id_field_name, sanitize, removeHs)
    elif input_format in ('smi', 'smiles'):
        suppl = __read_smiles(os.path.abspath(fname), sanitize)
    elif input_format == 'pkl':
        suppl = __read_pkl(os.path.abspath(fname))
    else:
        raise Exception("Input file format '%s' is not supported. It can be only sdf, sdf.gz, smi, pkl." % input_format)
    for mol, mol_name in suppl:
        yield mol, mol_name