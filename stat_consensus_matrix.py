import os
import pandas as pd


def filter_models(external_stat, filter=True):
    """
    select models for further analysis
    :param external_stat: model statistics
    :return: list of selected models
    """
    df = pd.read_csv(external_stat, sep='\t')
    if filter:
        df = df[(df['precision'] >= 0.8) | ((df['precision'] >= 0.7) & (df['EF'] > 3))]
        df = df[(df['TP'] > 2) & (df['num_uniq_F'] > 5) & (df['max_edge'] >= 8)]
    print('{} file: {} models'.format(os.path.basename(external_stat), df.shape[0]))
    return list(df['model'])


def build_matrix(path_screen, fmol, list_models):
    """
    Builds a binary matrix. Columns are model names, indexes are id compounds
    :param path_screen: path to files with screening results
    :param fmol: .smi file with active and inactive compounds used to create the database for virtual screening
    :param list_models: list of selected models
    :return: a binary matrix (dataframe) in which 1 - match between model and compounds,
                                                  0 - no match between model and compounds
    """
    nmols = list(pd.read_csv(fmol, sep='\t')['CID'])
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


def target_compounds_matrix(list_matrix):
    """
    Collects information about the number of matches between the compound and the target models
    :param list_matrix: list of dataframe binary matrices
    :return: a matrix on which columns are target names, indexes are id compounds
    """
    df_res = pd.DataFrame()
    for mmatrix in list_matrix:
        ntarget = mmatrix.index.name
        mmatrix.index.name = 'mol_name'
        mmatrix[ntarget] = mmatrix.sum(axis=1)

        df_res = pd.concat([df_res, mmatrix[[ntarget]]], axis=1)
    df_res = df_res.fillna(0)
    return df_res


def target_compounds_binary_matrix(df):
    df[df > 0] = 1
    return df


def main(pcompounds, pscreens, pstats, poutput, pouput_binary):

    if not os.path.exists(os.path.dirname(poutput)):
        os.mkdir(os.path.dirname(poutput))

    list_binary_matrices = []
    for ncom in os.listdir(pstats):
        nmark = ncom.split('.')[0].split('_')[-1]
        # find files of set
        pscreen = os.path.join(pscreens, nmark)
        ppstat = os.path.join(pstats, ncom)

        list_mod = filter_models(ppstat)
        df_t = build_matrix(pscreen, pcompounds, list_mod)
        list_binary_matrices.append(df_t)

    df_res = target_compounds_matrix(list_binary_matrices)
    df_res.to_csv(poutput, sep='\t')

    df = target_compounds_binary_matrix(df_res)
    df.to_csv(pouput_binary, sep='\t')


def calc_tanimoto(true_matrix, pred_matrix):
    df_res = pd.DataFrame(index=pred_matrix.index.tolist(), columns=['tanimoto', 'TP', 'FP'])
    for nmol in pred_matrix.index.tolist():
        try:
            df_res.at[nmol, 'tanimoto'] = sum(true_matrix.loc[nmol] & pred_matrix.loc[nmol]) / sum(
                true_matrix.loc[nmol] | pred_matrix.loc[nmol])
        except ZeroDivisionError:
            df_res.at[nmol, 'tanimoto'] = -1
        df_res.at[nmol, 'TP'] = sum(true_matrix.loc[nmol] & pred_matrix.loc[nmol])
        df_res.at[nmol, 'FP'] = sum(pred_matrix.loc[nmol]) - df_res.at[nmol, 'TP']
    return df_res


def starter_filter():
    main(pcompounds='data/compounds/assays10.act',
         pscreens='data/screen',
         pstats='data/models_stat',
         poutput='data/results/lessf4_filter_compounds_targets_matrix.csv',
         pouput_binary='data/results/lessf4_filter_compounds_targets_binary_matrix.csv')

    dfp = pd.read_csv('data/results/lessf4_filter_compounds_targets_binary_matrix.csv', sep='\t', index_col=0)
    dft = pd.read_csv('data/compounds/assays10.act', sep='\t', index_col=0)

    dfp = dfp.astype(int)

    cols = dfp.columns.tolist()
    dft = dft[cols]
    dfp = dfp[cols]

    df_res = calc_tanimoto(dft, dfp)
    df_res.to_csv('data/results/lessf4_filter_stat.txt', sep='\t')


def starter():
    main(pcompounds='data/compounds/assays10.act',
         pscreens='data/screen',
         pstats='data/models_stat',
         poutput='data/results/compounds_targets_matrix.csv',
         pouput_binary='data/results/compounds_targets_binary_matrix.csv')

    dfp = pd.read_csv('data/results/compounds_targets_binary_matrix.csv', sep='\t', index_col=0)
    dft = pd.read_csv('data/compounds/assays10.act', sep='\t', index_col=0)

    dfp = dfp.astype(int)

    cols = dfp.columns.tolist()
    dft = dft[cols]
    dfp = dfp[cols]

    df_res = calc_tanimoto(dft, dfp)
    df_res.to_csv('data/results/stat.txt', sep='\t')

starter_filter()
# starter()