import consensus_cdk2
import pandas as pd
import matplotlib.pyplot as plt

mol_act = 'test_set/active.smi'
mol_inact = 'test_set/inactive.smi'


def maker_path(status):
    a = ['screen/{}/2c6o'.format(status), 'pma/2c6o']
    b = ['screen/{}/2fvd'.format(status), 'pma/2fvd']
    c = ['screen/{}/2xmy'.format(status), 'pma/2xmy']
    d = ['screen/{}/5d1j'.format(status), 'pma/5d1j']
    return a, b, c, d


plog = 'CHEMBL301_class_log'
df_log = pd.read_csv(plog)
df_log = df_log[['cmp', 'pchembl_value', 'bioactivity_type']]
df_log = df_log[df_log['bioactivity_type'] == 'Ki']
df_log = df_log.rename(columns={'cmp': 'mol_name', 'pchembl_value': 'pKi'})
df_log = df_log.drop_duplicates(subset='mol_name')
df_log.index = list(df_log['mol_name'])
df_log.loc['num_features', :] = [0 for _ in range(df_log.shape[1])]

for n, pp in enumerate(maker_path('active')):
    pis = maker_path('inactive')
    tmol = pp[0].split('/')[-1]
    df_models_a = consensus_cdk2.pick_models(pp[1], pp[0])
    df_models_i = consensus_cdk2.pick_models(pp[1], pis[n][0])
    if df_models_a.empty and df_models_i.empty:
        continue

    df_models = df_models_a.append(df_models_i, sort=False)
    df_models = df_models.drop_duplicates('filename')
    df_cmod = pd.DataFrame(columns=list(df_models['filename']), index=['num_features'])
    df_cmod.loc['num_features', :] = list(df_models['count'])

    df_matrix_a = consensus_cdk2.maker_matrix(mol_act, df_models_a, 1)
    df_matrix_i = consensus_cdk2.maker_matrix(mol_inact, df_models_i, 0)

    df_matrix = pd.concat([df_matrix_a, df_matrix_i], ignore_index=False, sort=True)
    df_matrix = df_matrix.fillna(0)
    df_matrix = pd.merge(df_matrix, df_log[['pKi']], how='left', left_index=True, right_index=True)
    df_matrix = df_matrix.sort_values(by='y_pred', ascending=False)

    h = 'hashes/hashes_{}.txt'.format(tmol)
    df_hash = pd.read_csv(h, sep='\t')
    df_hash = df_hash[['filename', 'hash']]
    df_models_a = pd.merge(df_models_a, df_hash, on='filename', how='left')
    df_models_i = pd.merge(df_models_i, df_hash, on='filename', how='left')
    df_models_a = df_models_a.drop_duplicates(subset=['hash'])
    df_models_i = df_models_i.drop_duplicates(subset=['hash'])

    df_models = df_models_a.append(df_models_i, sort=False)
    df_models = df_models.drop_duplicates('filename')
    df_cmod_u = pd.DataFrame(columns=list(df_models['filename']), index=['num_features'])
    df_cmod_u.loc['num_features', :] = list(df_models['count'])

    df_matrix_a = consensus_cdk2.maker_matrix(mol_act, df_models_a, 1)
    df_matrix_i = consensus_cdk2.maker_matrix(mol_inact, df_models_i, 0)
    matrix_uni = pd.concat([df_matrix_a, df_matrix_i], ignore_index=False, sort=True)
    matrix_uni = matrix_uni.fillna(0)
    matrix_uni = pd.merge(matrix_uni, df_log[['pKi']], how='left', left_index=True, right_index=True)
    matrix_uni = matrix_uni.sort_values(by='y_pred', ascending=False)

    consensus_cdk2.maker_roc(tmol, df_matrix, 'consensus/roc_{}_full.png'.format(tmol))
    consensus_cdk2.maker_roc(tmol, matrix_uni, 'consensus/roc_{}_uniq.png'.format(tmol))

    df_matrix = pd.concat([df_cmod, df_matrix], ignore_index=False, sort=True)
    matrix_uni = pd.concat([df_cmod_u, matrix_uni], ignore_index=False, sort=True)
    df_matrix = df_matrix.sort_values(by='num_features', axis=1)
    matrix_uni = matrix_uni.sort_values(by='num_features', axis=1)
    df_matrix.index.name = 'mol_name'
    matrix_uni.index.name = 'mol_name'

    # for mm in df_models['filename']:
    #     for i in df_matrix[mm]:
    #         if (i - int(i)) != 0:
    #             print(mm, i)
    #
    # for mm in df_models['filename']:
    #     for i in matrix_uni[mm]:
    #         if (i - int(i)) != 0:
    #             print(mm, i)
    # for i in range(df_matrix.shape[0]):
    #     pred = df_matrix['y_pred'].iloc[i]
    #     try:
    #         if (pred - int(pred)) != 0:
    #             print(df_matrix.iloc[i])
    #     except ValueError:
    #         continue
    #
    # # for i in range(matrix_uni.shape[0]):
    # #     pred = matrix_uni['y_pred'].iloc[i]
    # #     try:
    # #         if (pred - int(pred)) != 0:
    # #             print(matrix_uni.iloc[i])
    # #     except ValueError:
    # #         continue
    #
    # for i in df_matrix.loc['num_features', :]:
    #     try:
    #         if (i - int(i)) != 0:
    #             print(i)
    #     except:
    #         print(i)
    # for i in matrix_uni.loc['num_features', :]:
    #     try:
    #         if (i - int(i)) != 0:
    #             print(i)
    #     except:
    #         print(i)

    df_matrix.to_csv('consensus/matrix_{}_full.csv'.format(tmol))
    matrix_uni.to_csv('consensus/matrix_{}_uniq.csv'.format(tmol))

    with open('consensus/{}_per_top.txt'.format(tmol), 'w') as f:
        f.write(tmol + '\n')
        f.write('per\tTP\tFP\tprecision\tEF\n')
        for ntop in [1, 2, 5, 10]:
            f.write(consensus_cdk2.calc_metrix(df_matrix, ntop))
        f.write('\nunique hashes\n')
        f.write('per\tTP\tFP\tprecision\tEF\n')
        for ntop in [1, 2, 5, 10]:
            f.write(consensus_cdk2.calc_metrix(matrix_uni, ntop))
