import os
import pandas as pd
import json
import sklearn.metrics as metrics
import matplotlib.pyplot as plt
plt.switch_backend('agg')
from collections import defaultdict


def pick_models(pp_pma, pp_screen):
    # dl = defaultdict(list)
    dl = {'filename': [], 'count': [], 'pp_screen': []}
    for i, pp in enumerate(os.listdir(pp_pma)):
        pma = os.path.join(pp_pma, pp)
        ll = os.path.join(pp_screen, pp.split('.')[0] + '.txt')
        with open(pma) as fpma:
            d = json.loads(fpma.readline().strip())
            num_uniq_features = len(set(tuple(feature[1]) for feature in d['feature_coords']))
            if (num_uniq_features >= 6) and (os.path.exists(ll)):
                dl['filename'].append(pp.split('.')[0].split('_')[1])
                dl['count'].append(num_uniq_features)
                dl['pp_screen'].append(ll)
    df = pd.DataFrame(dl)
    if df.empty:
        return df
    df = df.drop_duplicates(subset='filename')
    df = df.sort_values(by=['count'])
    return df


def pick_fit_models(pp_screen):
    models = []
    for pp in os.listdir(pp_screen):
        models.append(pp.split('.')[0])
    return models


def maker_matrix(pp_mols, df_models, status):
    if df_models.empty:
        return pd.DataFrame()
    df_mols = pd.read_csv(pp_mols, sep='\t')
    df_mols = df_mols.drop_duplicates(subset='mol_name')
    df = pd.DataFrame(columns=list(df_models['filename']), index=list(df_mols['mol_name']))
    for m, ll in zip(df_models['filename'], df_models['pp_screen']):
        with open(ll, 'r') as f:
            for mol in f.readlines():
                df.at[mol.strip(), m] = 1
    df = df.fillna(0)
    df.loc[:, 'y_pred'] = df.sum(axis=1)
    df['y_true'] = [status for _ in range(df.shape[0])]
    return df


def calc_metrix(df, const):
    df = df.sort_values(by=['y_pred'], ascending=False)
    n = round(df.shape[0]*(const/100))
    df_top = df[:n]
    TP = len(df_top[df_top['y_true'] == 1])
    FP = len(df_top[df_top['y_true'] == 0])
    all_active_mol = len(df[df['y_true'] == 1])
    all_inactive_mol = len(df[df['y_true'] == 0])
    precision = TP / (TP + FP)
    ef = (TP / (TP + FP)) / (all_active_mol / (all_inactive_mol + all_active_mol))
    return '{}\t{}\t{}\t{}\t{}\n'.format(const, TP, FP, precision, ef)


def maker_roc(title, df_res, path_out):
    df_res = df_res.sort_values(by=['y_pred'], ascending=False)
    fpr, tpr, threshold = metrics.roc_curve(list(df_res['y_true']), list(df_res['y_pred']))
    roc_auc = metrics.auc(fpr, tpr)
    plt.title(title)
    plt.plot(fpr, tpr, 'b', label='AUC = %0.2f' % roc_auc)
    plt.legend(loc='lower right')
    plt.plot([0, 1], [0, 1], 'r--')
    plt.xlim([0, 1])
    plt.ylim([0, 1])
    plt.ylabel('True Positive Rate')
    plt.xlabel('False Positive Rate')
    plt.savefig(path_out)
    plt.close()
