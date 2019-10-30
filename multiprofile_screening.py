import os
import numpy as np
import pandas as pd
from screen_db import screen_db, get_comp_names_from_db


def build_matrix(path_screen, nmols, list_models):
    """
    Builds a binary matrix. Columns are model names, indexes are id compounds
    :param path_screen: path to files with screening results
    :param fmol: .smi file with active and inactive compounds used to create the database for virtual screening
    :param list_models: list of selected models
    :return: a binary matrix (dataframe) in which 1 - match between model and compounds,
                                                  0 - no match between model and compounds
    """
    df = pd.DataFrame(columns=list_models, index=nmols)
    ntarget = os.path.split(path_screen)[1]
    df.index.name = ntarget
    df = df.fillna(0)
    for num_db in os.listdir(path_screen):
        for ff in os.listdir(os.path.join(path_screen, num_db)):
            if ff.split('.')[0] in list_models:
                pscreenfile = os.path.join(path_screen, num_db, ff)
                match_compounds = [int(mol.strip()) for mol in open(pscreenfile).readlines()]
                for compound in match_compounds:
                    df.at[compound, ff.split('.')[0]] = 1
    df = df.fillna(0)
    return df


def multiprofile_screening(path_db, target_csv, path_screens, ncpu):
    path_chembls = '/home/akutlushina/chembl'
    list_targets = [ii.strip().split('\t')[0] for ii in open(target_csv).readlines() if ii[0] != '#']

    for chembl_id in list_targets:
        path_screen = os.path.join(path_screens, chembl_id)
        path_query = os.path.join(path_chembls, chembl_id, 'models')
        screen_db(db_fname=path_db,
                  queries=path_query,
                  output=path_screen,
                  input_sdf=None,
                  match_first_conf=True,
                  ncpu=ncpu)

        nmols = get_comp_names_from_db(path_db)
        list_models = os.listdir(path_query)
        df_matrix = build_matrix(path_screen, nmols, list_models)
        df_matrix.to_csv(os.path.join(path_screen, '{}_matrix.csv'.format(chembl_id)))

        df_matrix.loc[:, 'sum'] = df_matrix.sum(axis=1)
        df_matrix['status'] = np.where(df_matrix['sum'] > 0, 1, 0)
        df_matrix[['status']].to_csv(os.path.join(path_screen, '{}_predicted.csv'.format(chembl_id)))
