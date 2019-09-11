import os

import pandas as pd
from rdkit.Chem import Conformer, Descriptors3D
from rdkit.Geometry import Point3D
from rdkit import Chem
from pharmacophore import Pharmacophore
from collections import defaultdict


def calc_descr(pmol):
    return Descriptors3D.Asphericity(pmol), Descriptors3D.Eccentricity(pmol), \
           Descriptors3D.InertialShapeFactor(pmol), Descriptors3D.NPR1(pmol), \
           Descriptors3D.NPR2(pmol), Descriptors3D.PMI1(pmol), \
           Descriptors3D.PMI2(pmol), Descriptors3D.PMI3(pmol),\
           Descriptors3D.RadiusOfGyration(pmol), Descriptors3D.SpherocityIndex(pmol)


def add_descr():
    df_all = pd.DataFrame()
    pp = '/home/akutlushina/compounds/chembl_test'

    for mdls in os.listdir(pp):
        ppath = os.path.join(pp, mdls, 'binding')

        df = pd.DataFrame()

        for fls in os.listdir(ppath):
            if fls.split('_')[0] == 'external':
                df_one = pd.read_csv(os.path.join(ppath, fls), sep='\t')
                df = df.append(df_one)

        if not df.empty:
            df.to_csv('/home/akutlushina/compounds/res/external_stat_{}_{}.txt'.format(mdls, 'binding'), index=None, sep='\t')
            df['ID'] = [mdls for _ in range(df.shape[0])]
            df['categories'] = ['binding' for _ in range(df.shape[0])]
            df_all = df_all.append(df)
    return df_all


def collected():
    pp = os.path.abspath('/home/akutlushina/compounds/res/')
    out = open('/home/akutlushina/compounds/external_stat.csv', 'w')

    for fs in os.listdir(pp):
        out.write(fs)
        with open(os.path.join(pp, fs), 'r') as f:
            for line in f.readlines():
                out.write(line)
            out.write('\n')


def max_edge(coords):
    edge = 0
    for i, c1 in enumerate(coords):
        for j, c2 in enumerate(coords[i+1:]):
            e = ((c1[1][0] - c2[1][0])**2 + (c1[1][1] - c2[1][1])**2 + (c1[1][2] - c2[1][2])**2) ** (1/2)
            if e > edge:
                edge = e    
    return e


def mol_h(query):
    q = Pharmacophore(cached=True)
    q.load_from_pma(query)

    pmol = Chem.RWMol()
    all_coords = q.get_feature_coords(ids=None)
    edge = max_edge(all_coords)
    for _ in all_coords:
        a = Chem.Atom(1)
        pmol.AddAtom(a)
    c = Conformer(len(all_coords))
    for i, coords in enumerate(all_coords):
        c.SetAtomPosition(i, Point3D(*coords[1]))
    pmol.AddConformer(c, True)
    return pmol, edge


pp = '/home/akutlushina/compounds/chembl_test/'

d = defaultdict(list)
for mdls in os.listdir(pp):
    ppath = os.path.join(pp, mdls, 'binding', 'models')

    if os.path.exists(ppath):
        for fls in os.listdir(ppath):
            for pma in os.listdir(os.path.join(ppath, fls)):

                pmol, edge = mol_h(os.path.join(ppath, fls, pma))
                Asphericity, Eccentricity, InertialShapeFactor, \
                NPR1, NPR2, PMI1, PMI2, PMI3, RadiusOfGyration, SpherocityIndex = calc_descr(pmol)
                d['model'].append(pma)
                d['ID'].append(mdls)
                d['categories'].append('binding')
                d['max_edge'].append(edge)
                d['Asphericity'].append(Asphericity)
                d['Eccentricity'].append(Eccentricity)
                d['InertialShapeFactor'].append(InertialShapeFactor)
                d['NPR1'].append(NPR1)
                d['NPR2'].append(NPR2)
                d['PMI1'].append(PMI1)
                d['PMI2'].append(PMI2)
                d['PMI3'].append(PMI3)
                d['RadiusOfGyration'].append(RadiusOfGyration)
                d['SpherocityIndex'].append(SpherocityIndex)

df = pd.DataFrame(d)
if not df.empty:
    df.to_csv('/home/akutlushina/compounds/res/model_stat_binding.txt', index=None, sep='\t')

df = add_descr()
df.to_csv('/home/akutlushina/compounds/external_stat_binding.txt', index=None, sep='\t')
