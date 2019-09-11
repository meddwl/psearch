#!/usr/bin/env python3

import os
import argparse
from rdkit import Chem
from rdkit.Chem import ChemicalFeatures

import xyz_pp
from pharmacophore import Pharmacophore


def main(path_set, conf_sdf, out_comp, out_model, process_factory):
    query_path = os.path.join(path_set, 'models')
    path_screen = os.path.join(path_set, 'screen')

    for m in os.listdir(query_path):
        for model in os.listdir(os.path.join(query_path, m)):
            model_pma = os.path.join(query_path, m, model)
            screen_act = os.path.join(path_screen, m, 'screen_active_{}.txt'.format(model.split('.')[0]))

            q = Pharmacophore(cached=True)
            q.load_from_pma(model_pma)

            pmol = q.get_mol()
            pmol.SetProp('_Name', model.split('.')[0])
            writer = Chem.SDWriter(os.path.join(out_model, '{}.sdf'.format(model.split('.')[0])))
            writer.write(pmol)

            supp = xyz_pp.main(query_fname=model_pma,
                               mol_sdf=conf_sdf,
                               screen_model=screen_act,
                               process_factory=process_factory,
                               db_bin_step=1,
                               get_transform_matrix=True)
            xyz_pp.save_sdf(mols=supp, out=os.path.join(out_comp, '{}_compounds.sdf'.format(model.split('.')[0])))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Screen SQLite DB with compounds against pharmacophore queries. '
                                                 'Output is printed out in STDOUT.',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-i', '--path_set', required=True,
                        help='')
    parser.add_argument('-conf', '--conformers', metavar='conformers.sdf',
                        help='SDF file with conformers structures.')
    parser.add_argument('-oc', '--output_compounds', metavar='output.sdf',
                        help='output sdf file with transform matrix of molecule.')
    parser.add_argument('-om', '--output_model', metavar='model.sdf',
                        help='output sdf file with model.')
    parser.add_argument('--fdef', default=os.path.join(os.path.split(os.getcwd())[0], 'pmapper/smarts_features.fdef'),
                        help='fdef-file with pharmacophore feature definition.')

    args = vars(parser.parse_args())
    for o, v in args.items():
        if o == "path_set": path_set = v
        if o == "conformers": conf_sdf = v
        if o == "output_compounds": out_comp = v
        if o == "output_model": out_model = v
        if o == "fdef": fdef_fname = v

    process_factory = ChemicalFeatures.BuildFeatureFactory(fdef_fname) if fdef_fname else None
    for outd in [out_comp, out_model]:
        if not os.path.exists(outd):
            os.mkdir(outd)

    main(path_set=path_set,
         conf_sdf=conf_sdf,
         out_comp=out_comp,
         out_model=out_model,
         process_factory=process_factory)
