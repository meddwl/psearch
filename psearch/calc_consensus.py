import os
import argparse


def cal_consensus(mols_screen, act_mol, inact_mol):
    tp = []
    fp = []
    for ll in mols_screen:
        tp = tp + [ii.split('\t')[1] for ii in open(ll['active']).readlines()]
        if os.path.exists(ll['inactive']):
            fp = fp + [ii.split('\t')[1] for ii in open(ll['inactive']).readlines()]

    p = len(open(act_mol).readlines())
    n = len(open(inact_mol).readlines())
    fn = p - tp
    tn = n - fp
    precision = round(tp / (tp + fp), 3)
    recall = round(tp / (tp + fn), 3)
    fpr = round(fp / (tn + fp), 3)
    f1 = round((2 * precision * recall) / (precision + recall), 3)
    f2 = round((5 * precision * recall) / (4 * precision + recall), 3)
    f05 = round((1.25 * precision * recall) / (0.25 * precision + recall), 3)
    ef = round((tp / (tp + fp)) / (p / (p + n)), 3)

    return len(set(tp)), len(set(fp)), precision, recall, fpr, f1, f2, f05, ef



def main(act_mol, inact_mol, path_screen, out_path):
    s1_models = []
    s2_models = []
    for ff in os.listdir(os.path.join(path_screen, 'active')):
        if ff.split('_')[0] == 'centroid':
            s1_models.append({'active': os.path.join(os.path.abspath(path_screen), 'active', ff),
                                     'inactive': os.path.join(os.path.abspath(path_screen), 'inactive', ff)})
        else:
            s2_models.append({'active': os.path.join(os.path.abspath(path_screen), 'active', ff),
                                     'inactive':os.path.join(os.path.abspath(path_screen), 'inactive', ff)})

    with open(out_path, 'w') as f:
        f.write('strategy 1\t' + '\t'.join(map(str, cal_consensus(s1_models, act_mol, inact_mol))) + '\n')
        f.write('strategy 2\t' + '\t'.join(map(str, cal_consensus(s2_models, act_mol, inact_mol))) + '\n')
        f.write('all models\t' + '\t'.join(map(str, cal_consensus(s1_models + s2_models, act_mol, inact_mol))) + '\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='select active and inactive compounds'
                                                 'based on given values (act_threshold and inact_threshold)')
    parser.add_argument('-s', '--path_screen', metavar='input.smi', required=True,
                        help='input SMILES file name. It should contain three columns separated by whitespaces: '
                             'SMILES, name, activity. No header.')
    parser.add_argument('-am', '--active_molecules', metavar='', required=True,
                        help='')
    parser.add_argument('-im', '--inactive_molecules', metavar='', required=True,
                        help='')
    parser.add_argument('-o', '--output_path', metavar='', default=None,
                        help='')

    args = vars(parser.parse_args())
    for o, v in args.items():
        if o == "path_screen": path_screen = v
        if o == "active_molecules": act_mol = v
        if o == "inactive_molecules": inact_mol = v
        if o == "output_path": output_path = v

    if not output_path:
        output_path = os.path.join(os.path.dirname(act_mol), 'consensus_results.txt')
    if not os.path.exists(os.path.dirname(output_path)):
        os.makedirs(output_path)

    main(path_screen=path_screen,
         act_mol=act_mol,
         inact_mol=inact_mol,
         out_path=output_path)