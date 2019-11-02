#!/usr/bin/env python3
# author          : Pavel Polishchuk
# date            : 13.07.16
# license         : BSD-3
#==============================================================================

import os
import sys
import marshal
import argparse
import sqlite3 as lite

from rdkit import Chem
from multiprocessing import cpu_count, Pool
from .read_input import read_input
from pmapper.pharmacophore import Pharmacophore
from pmapper.customize import load_smarts, load_factory
from pmapper.utils import load_multi_conf_mol


def create_tables(cursor, bin_step, smarts):
    cursor.execute("CREATE TABLE conformers(conf_id INTEGER PRIMARY KEY AUTOINCREMENT, mol_name TEXT, "
                   "mol_stereo_id TEXT, fp BLOB)")
    cursor.execute("CREATE INDEX mol_name_idx ON conformers(mol_name)")
    cursor.execute("CREATE TABLE settings(bin_step NUMERIC)")
    cursor.execute("INSERT INTO settings VALUES(?)", (bin_step, ))
    cursor.execute("CREATE TABLE smarts_features(label TEXT, value TEXT)")
    if smarts:
        for k, values in smarts.items():
            for v in values:
                cursor.execute("INSERT INTO smarts_features VALUES(?, ?)", (k, Chem.MolToSmarts(v)))
    cursor.execute("CREATE TABLE feature_coords(conf_id INTEGER, feature_label TEXT, "
                   "x NUMERIC, y NUMERIC, z NUMERIC, "
                   "FOREIGN KEY(conf_id) REFERENCES conformers(conf_id))")


def insert_res_db(cur, res, stereo_id):
    for item in res:   # mol_name, coords, fp
        if stereo_id:
            mol_name, mol_stereo_id = item[0].rsplit("_", 1)
        else:
            mol_name = item[0]
            mol_stereo_id = '0'
        data = [mol_name, mol_stereo_id, item[2]]
        cur.execute("INSERT INTO conformers VALUES(NULL, ?, ?, ?)", data)
        # store coords
        cur.execute("SELECT MAX(conf_id) FROM conformers")
        conf_id = cur.fetchone()[0]
        cur.executemany("INSERT INTO feature_coords VALUES (?, ?, ?, ?, ?)", ((conf_id, i[0], *i[1]) for i in item[1]))


def compress_db(cursor, store_coords):
    # remove duplicated hashes for the same mol
    cursor.execute("DELETE FROM conformers WHERE rowid NOT IN (SELECT MIN(rowid) "
                   "FROM conformers GROUP BY mol_name, mol_stereo_id, ph_hash)")
    if store_coords:
        cursor.execute("DELETE FROM feature_coords WHERE conf_id NOT IN (SELECT conf_id "
                       "FROM conformers)")


def insert_res_txt(f, res, lines_set, stereo_id):
    # function is called if hash is not None, no need for additional check
    for item in res:   # mol_name, hash, coords, fp
        if stereo_id:
            mol_name, mol_stereo_id = item[0].rsplit("_", 1)
        else:
            mol_name = item[0]
            mol_stereo_id = '0'
        record = (mol_name, mol_stereo_id, item[1])
        if record not in lines_set:
            f.write("%s\t%s\t%s\n" % record)
            lines_set.add(record)


def prep_input(fname, id_field_name, smarts, bin_step, multiconf):
    if fname is None:
        read_iterator = read_input(fname, 'sdf', id_field_name, removeHs=False)
    else:
        read_iterator = read_input(fname, id_field_name=id_field_name, removeHs=False)
    for mol, mol_name in read_iterator:
        yield mol, mol_name, smarts, bin_step, multiconf


def map_process_mol(args):
    return process_mol(*args)


def process_mol(mol, mol_name, smarts, bin_step, multiconf):
    # process_factory is within process scope only
    if multiconf:
        if 'process_factory' in globals():
            ps = load_multi_conf_mol(mol, factory=process_factory, bin_step=bin_step)
        else:
            ps = load_multi_conf_mol(mol, smarts_features=smarts, bin_step=bin_step)
        output = []
        for p in ps:
            coords = p.get_feature_coords()
            fp_bin = marshal.dumps(p.get_fp())
            output.append((mol_name, coords, fp_bin))
        return output
    else:
        p = Pharmacophore(bin_step, cached=True)
        if 'process_factory' in globals():
            p.load_from_feature_factory(mol, process_factory)
        elif smarts:
            p.load_from_smarts(mol, smarts)
        coords = p.get_feature_coords()
        fp_bin = marshal.dumps(p.get_fp())
        return [(mol_name, coords, fp_bin)]


def pool_init(fdef_fname):
    global process_factory
    process_factory = load_factory(fdef_fname)


def main_params(conformers_fname, dbout_fname, bin_step, rewrite_db, id_field_name, stereo_id,
                smarts_features_fname, rdkit_factory, ncpu, verbose):

    # check DB existence
    if dbout_fname is not None and os.path.isfile(dbout_fname):
        if rewrite_db:
            os.remove(dbout_fname)
        else:
            raise FileExistsError("DB exists. To rewrite it add -r key to the command line call or remove "
                                  "database manually.")

    # load smarts features
    if rdkit_factory is None:
        if smarts_features_fname is not None:
            smarts = load_smarts(smarts_features_fname)
        else:
            smarts = load_smarts()
    else:
        smarts = None

    multiconf = conformers_fname.lower().endswith('.pkl')

    # open db
    conn = lite.connect(dbout_fname)
    cur = conn.cursor()
    create_tables(cur, bin_step, smarts)

    nprocess = max(min(ncpu, cpu_count()), 1)

    if smarts:
        p = Pool(nprocess)
    else:
        p = Pool(nprocess, initializer=pool_init, initargs=[rdkit_factory])

    lines_set = set()

    counter = 0

    try:
        for i, res in enumerate(p.imap_unordered(map_process_mol,
                                                 prep_input(conformers_fname, id_field_name, smarts, bin_step,
                                                            multiconf),
                                                 chunksize=10), 1):
            counter += len(res)
            if dbout_fname is not None:
                insert_res_db(cur, res, stereo_id)
                if i % 1000 == 0:
                    conn.commit()
            if verbose and i % 100 == 0:
                sys.stderr.write('\rprocessed %i molecules and %i conformers/pharmacophores' % (i, counter))
                sys.stderr.flush()

        if dbout_fname is not None:
            # compress_db(cur, store_coords)
            conn.commit()

    finally:
        p.close()
        if dbout_fname is not None:
            conn.close()

    if verbose:
        sys.stderr.write("\n")


def entry_point():
    parser = argparse.ArgumentParser(description='Create DB with pharmacophore representation of input compound '
                                                 'library.', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--input', metavar='conformers.sdf', required=False, default=None,
                        help='SDF or SDF.GZ file with conformers structures. '
                             'Conformers of the same molecule should have identical names. '
                             'Alternatively PKL files can be used as input, which are stored pickled tuples of '
                             'multi-conformer molecules and their names. '
                             'If omitted STDIN will be parsed as SDF file.')
    parser.add_argument('-d', '--dbout', metavar='output.db', required=True,
                        help='output DB SQLite file.')
    parser.add_argument('-b', '--bin_step', default=1,
                        help='binning step. Default: 1.')
    parser.add_argument('-s', '--smarts_features', metavar='smarts.txt', default=None,
                        help='text file with definition of pharmacophore features. Tab-separated two columns: '
                             'first - SMARTS, second - feature name. If file name is not specified the default SMARTS '
                             'patterns will be used.')
    parser.add_argument('--rdkit_factory', metavar='features.fdef', default=None,
                        help='text file with definition of pharmacophore features in RDKit format. If file name is not '
                             'specified the default SMARTS patterns will be used. This option has a priority over '
                             'smarts_features if both were specified.')
    parser.add_argument('-f', '--id_field_name', metavar='field_name', default=None,
                        help='field name of compound ID (sdf). If omitted for sdf molecule titles will be used or '
                             'auto-generated names. Please note if you use SDTIN as input, molecule names should be '
                             'stored in a title field of a MOL block, property fields will not be read.')
    parser.add_argument('--stereo_id', action='store_true', default=False,
                        help='set this option if mol names contain stereo_id after last underscore character '
                             '(e.g. MolName_1, Mol_name_2, etc). Then different stereoisomers will be considered '
                             'together when duplicated pharmacophores will be removed. Otherwise each stereoisomer '
                             'will be treated as an individual compound.')
    parser.add_argument('-r', '--rewrite', action='store_true', default=False,
                        help='rewrite existed DB.')
    parser.add_argument('-c', '--ncpu', default=1,
                        help='number of cpu used. Default: 1.')
    parser.add_argument('-v', '--verbose', action='store_true', default=False,
                        help='print progress.')

    args = vars(parser.parse_args())
    for o, v in args.items():
        if o == "input": conformers_fname = v
        if o == "dbout": dbout_fname = v
        if o == "bin_step": bin_step = float(v)
        if o == "smarts_features": smarts_features = v
        if o == "id_field_name": id_field_name = v
        if o == "rewrite": rewrite_db = v
        if o == "verbose": verbose = v
        if o == "ncpu": ncpu = int(v)
        if o == "stereo_id": stereo_id = v
        if o == "rdkit_factory": rdkit_factory = v

    main_params(dbout_fname=dbout_fname,
                smarts_features_fname=smarts_features,
                rdkit_factory=rdkit_factory,
                conformers_fname=conformers_fname,
                bin_step=bin_step,
                rewrite_db=rewrite_db,
                id_field_name=id_field_name,
                stereo_id=stereo_id,
                verbose=verbose,
                ncpu=ncpu)


if __name__ == '__main__':
    entry_point()
