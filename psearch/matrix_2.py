import argparse
import os
import pandas as pd


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


def main(path_screen, path_act, path_inact, path_all_models, act_duplicate, inact_duplicate, out_path):
    acts = [i.strip().split('\t')[1] for i in open(path_act).readlines()]
    inacts = [i.strip().split('\t')[1] for i in open(path_inact).readlines()]

    if act_duplicate:
        act_duplicate = [i.strip().split('\t')[1] for i in open(act_duplicate).readlines()]
        inact_duplicate = [i.strip().split('\t')[1] for i in open(inact_duplicate).readlines()]

        acts = list(set(acts).difference(act_duplicate))
        inacts = list(set(inacts).difference(inact_duplicate))

    nmols = acts + inacts

    list_models = [i.split('.')[0] for i in os.listdir(path_all_models)]
    df = build_matrix(path_screen, nmols, list_models)
    df.to_csv(out_path)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='select active and inactive compounds'
                                                 'based on given values (act_threshold and inact_threshold)')
    parser.add_argument('-s', '--path_screen', metavar='input.smi', required=True,
                        help='')
    parser.add_argument('-m', '--path_all_models', metavar='input.smi', required=True,
                        help='')
    parser.add_argument('-pa', '--path_act', metavar='', required=True,
                        help='')
    parser.add_argument('-pi', '--path_inact', metavar='', required=True,
                        help='')
    parser.add_argument('-ad', '--act_duplicate', metavar='active.smi', default=None,
                        help='')
    parser.add_argument('-id', '--inact_duplicate', metavar='inactive.smi', default=None,
                        help='')
    parser.add_argument('-o', '--output_path', metavar='', default=None,
                        help='')

    args = vars(parser.parse_args())
    for o, v in args.items():
        if o == "path_screen": path_screen = v
        if o == "path_all_models": path_all_models = v
        if o == "path_act": path_act = v
        if o == "path_inact": path_inact = v
        if o == "act_duplicate": act_duplicate = v
        if o == "inact_duplicate": inact_duplicate = v
        if o == "output_path": output_path = v

    if act_duplicate:
        output_path = os.path.join(os.path.dirname(os.path.abspath(path_screen)),
                                   'matrix/{}_{}'.format(act_duplicate.split('_')[-2], act_duplicate.split('_')[-1]))
    else:
        output_path = os.path.join(os.path.dirname(os.path.abspath(path_screen)),
                                   'matrix/{}.txt'.format(os.path.abspath(path_screen).split('/')[-2])

    if not os.path.exists(os.path.dirname(output_path)):
        os.makedirs(output_path)

    main(path_screen=path_screen,
         path_act=path_act,
         path_inact=path_inact,
         path_all_models=path_all_models,
         act_duplicate=act_duplicate,
         inact_duplicate=inact_duplicate,
         out_path=output_path)