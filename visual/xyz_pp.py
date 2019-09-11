#!/usr/bin/env python3

import os
import sys
import argparse
import pandas as pd
from collections import defaultdict
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import ChemicalFeatures

from read_input import read_input
from pharmacophore import Pharmacophore


def compare_fp(query_fp, fp):
    return (query_fp & fp) == query_fp

def main(query_fname, mol_sdf, screen_model, process_factory, db_bin_step, get_transform_matrix=True):
    # load query
    q = Pharmacophore()
    q.load_from_pma(query_fname)
    q_fp = q.get_fp()

    # check bin steps
    if q.get_bin_step() != db_bin_step:
        sys.stderr.write('Model has a different bin step from compounds in database. It would be skipped.\n')
        raise Exception('Incompatible bin step')

    #suppl = Chem.SDMolSupplier(mol_sdf)
    if os.path.getsize(screen_model) == 0:
        return '{} is empty'.format(screen_model)
    df_screen = pd.read_csv(screen_model, sep='\t', header=None)
    mol_names = list(set(df_screen[0]))

    mol_dict = defaultdict(list)
    read_iterator = read_input(mol_sdf, id_field_name=None, removeHs=False)
    mname = ''
    for mol, mname in read_iterator:
        mname = mname.split('_')[0]
        if mname in mol_names:
            mol_dict[mname].append(mol)
    print(os.path.basename(query_fname), len(mol_dict[mname]))
    suppl_out = []
    for mol_name, mols in mol_dict.items():
        for mol in mols:
            db_p = Pharmacophore(db_bin_step)
            db_p.load_from_feature_factory(mol, process_factory)
            feature_coords = db_p.get_feature_coords()

            p = Pharmacophore()
            p.load_from_feature_coords(feature_coords)
            res = p.fit_model(q, get_transform_matrix=get_transform_matrix)

            if res:
                AllChem.TransformMol(mol, res[1])
                mol.SetProp('transformation_matrix', ', '.join(str(i) for i in res[1].reshape((16,))))
                suppl_out.append(mol)
                res_string = '\t'.join(map(str, (mol_name, ','.join(map(str, res[0]))))) + '\n'
                sys.stdout.write(res_string)
                break #catch one conformer of molecule
    return suppl_out

def save_sdf(mols, out):
    if type(mols) is str:
        print(mols, '\n')
    else:
        writer = Chem.SDWriter(out)
        for mol in mols:
            writer.write(mol)
        writer.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Screen SQLite DB with compounds against pharmacophore queries. Output is printed out in STDOUT.',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-q', '--query_fname', metavar='model.pma', required=True,
                        help='pharmacophore models. Several models can be specified. Model format is recognized by file extension. Currently pma (pmapper) and pml (LigandScout) formats are supported.')
    parser.add_argument('-m', '--mol_sdf', metavar='molecules.sdf', required=True,
                        help='')
    parser.add_argument('-s', '--screen_model', metavar='screen_mol.txt', required=True,
                        help='path to output text file with screening results.')
    parser.add_argument('-tm', '--get_transform_matrix', default=True,
                        help='')
    parser.add_argument('-o', '--output', metavar='output.sdf', required=True,
                        help='output sdf file with transform matrix of molecule.')
    parser.add_argument('-b', '--bin_step', default=1,
                        help='binning step. Default: 1.')
    parser.add_argument('-f', '--fdef_fname', default=os.path.join(os.path.split(os.getcwd())[0], 'pmapper/smarts_features.fdef'),
                        help='fdef-file with pharmacophore feature definition.')

    args = vars(parser.parse_args())
    for o, v in args.items():
        if o == "query_fname": query_fname = v
        if o == "mol_sdf": mol_sdf = v
        if o == "screen_model": screen_model = v
        if o == "get_transform_matrix": get_transform_matrix = v
        if o == "output": out = v
        if o == "bin_step": bin_step = float(v)
        if o == "fdef_fname": fdef_fname = v

    process_factory = ChemicalFeatures.BuildFeatureFactory(fdef_fname) if fdef_fname else None

    suppl_out = main(query_fname=query_fname,
                     mol_sdf=mol_sdf,
                     screen_model=screen_model,
                     get_transform_matrix=get_transform_matrix,
                     db_bin_step=bin_step,
                     process_factory=process_factory)

    if os.path.exists(out):
        save_sdf(suppl_out, out)



