#!/usr/bin/env python3
# author          : Pavel Polishchuk
# date            : 12.05.16
# license         : BSD-3
# ==============================================================================

import os
import sys
import argparse
import marshal
import sqlite3 as lite
from collections import defaultdict
from pharmacophore import Pharmacophore


def compare_fp(query_fp, fp):
    return (query_fp & fp) == query_fp


def load_filtered_confs(cur, q_fp):
    cur.execute("SELECT DISTINCT(mol_name) FROM conformers")
    mol_names = [i[0] for i in cur.fetchall()]
    for mol_name in mol_names:
        # load fp for all conformers
        cur.execute("SELECT conf_id, fp FROM conformers WHERE mol_name = ?", (mol_name,))
        data = cur.fetchall()

        conf_ids = []
        for conf_id, fp in data:
            if compare_fp(q_fp, marshal.loads(fp)):
                conf_ids.append(conf_id)
        yield mol_name, conf_ids


def load_pharmacophores(cur, conf_ids):
    # load pharmacophores for selected conformers
    sql = "SELECT conf_id, feature_label, x, y, z FROM feature_coords WHERE conf_id IN (%s)" % \
          ','.join(['?'] * len(conf_ids))
    cur.execute(sql, conf_ids)
    res = cur.fetchall()
    confs = defaultdict(list)
    for r in res:
        confs[r[0]].append((r[1], tuple(r[2:])))
    return confs


def main(dbs_fname, path_pma, path_screen):
    if path_screen is None:
        path_screen = os.path.join(os.path.split(path_pma)[0], 'screen')

    for query_fname in os.listdir(path_pma):
        for in_db in dbs_fname:
            out_fname = os.path.join(path_screen, 'screen_{}_{}.txt'.format(
                os.path.basename(in_db).split('.')[0], query_fname.split('.')[0]))
            
            with open(out_fname, 'w') as out_f:
                q = Pharmacophore(cached=True)
                q.load_from_pma(os.path.join(path_pma, query_fname))
                q_fp = q.get_fp()
    
                conn = lite.connect(in_db)
                cur = conn.cursor()
    
                cur.execute("SELECT bin_step FROM settings")
                db_bin_step = cur.fetchone()[0]
                # check bin steps
                if q.get_bin_step() != db_bin_step:
                    sys.stderr.write('Model has a different bin step from compounds in database. '
                                     'It would be skipped.\n')
                    raise Exception('Incompatible bin step')

                for mol_name, conf_ids in load_filtered_confs(cur, q_fp):
                    if conf_ids:
                        for conf_id, feature_coords in load_pharmacophores(cur, conf_ids).items():
                            p = Pharmacophore(bin_step=db_bin_step, cached=True)
                            p.load_from_feature_coords(feature_coords)
                            res = p.fit_model(q)
                            if res:
                                res_string = '\t'.join(map(str, (mol_name, conf_id, ','.join(map(str, res))))) + '\n'
                                out_f.write(res_string)
                                break  # match the first conformer


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Screen SQLite DB with compounds against pharmacophore queries. '
                                                 'Output is printed out in STDOUT.',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-d', '--databases', metavar='input_active.db, input_inactive.db',
                        nargs='+', type=int, required=True,
                        help='SQLite input files (example active and inactive) with compounds to screen.')
    parser.add_argument('-q', '--query_path', metavar='model.pma', required=True,
                        help='pharmacophore models. Several models can be specified. Model format is recognized by '
                             'file extension. Currently pma (pmapper) and pml (LigandScout) formats are supported.')
    parser.add_argument('-o', '--output_path', default=None,
                        help='path to output text file with screening results.')

    args = vars(parser.parse_args())
    for o, v in args.items():
        if o == "databases": dbs_fname = v
        if o == "query_path": query_path = v
        if o == "output_path": out_fname = v

    main(dbs_fname=dbs_fname,
         path_pma=query_path,
         path_screen=out_fname)
