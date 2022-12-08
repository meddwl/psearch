#!/usr/bin/env python3
import os
import sys
import argparse

from math import sqrt
from itertools import combinations

from pmapper.pharmacophore import Pharmacophore as P


def get_max_dist(p):
    coords = p.get_feature_coords()
    dist = 0
    for a, b in combinations(coords, 2):
        d = sqrt((a[1][0] - b[1][0]) ** 2 + (a[1][1] - b[1][1]) ** 2 + (a[1][2] - b[1][2]) ** 2)
        if d > dist:
            dist = d
    return round(dist, 2)


def get_num_features(p):
    coords = p.get_feature_coords()
    return len(coords)


def get_num_distinct_features(p):
    coords = [i[1] for i in p.get_feature_coords()]
    return len(set(coords))


def get_num_polarfeature(p):
    nf = p.get_features_count()
    n = nf['A'] + nf['D'] + nf['P'] + nf['N']
    return n


def get_labels(p):
    coords = p.get_feature_coords()
    l = ''.join([i[0] for i in coords])
    return l


def calc_pharm_graph_desc(p):
    labels = get_labels(p)
    nf_dist = get_num_distinct_features(p)
    nf_polar = get_num_polarfeature(p)
    nf = get_num_features(p)
    dist = get_max_dist(p)
    return dist, nf, nf_dist, nf_polar, labels


def create_parser():
    parser = argparse.ArgumentParser(description='Get information about a pharmacophore graphs:\n'
                                                 'pharm_id - pharmacophore id;\t'
                                                 'max_dist - maximum distance between features in pharmacophore graph;\n'
                                                 'nf - number of features in pharmacophore graph;\n'
                                                 'nf_dist - number of distinct features in pharmacophore graph;\n'
                                                 'nf_polar - number of polar features (N, P, A, D) in pharmacophore graph;\n'
                                                 'labels - listed features of pharmacophore graph, \t'
                                                 'where a: aromatic centers, D: hydrogen donor centers, '
                                                 'A: hydrogen acceptor centers, P: positive charged centers, '
                                                 'N: negative charged centers, H: hydrophobic centers',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--models', type=str, required=True,
                        help='Path to a folder with pharmacophore models (only .pma or .xyz formats) or '
                             'to a pharmacophore model file.')
    parser.add_argument('-o', '--output', metavar='FILENAME', required=False, default=None,
                        help='Output path to a .txt file where will be saved pharmacophore graphs properties.'
                             'If omitted, the path will be generated automatically relative to project directory.')
    return parser


def entry_point():
    parser = create_parser()
    args = parser.parse_args()

    output = args.output if args.output else os.path.join(os.path.dirname(args.models), 'pharmacophore_properties.txt')

    if os.path.isdir(args.models):
        list_models = [os.path.join(args.models, i) for i in os.listdir(args.models)]
    else:
        list_models = [args.models]

    w = open(output, 'w')
    w.write('pharm_id\tmax_dist\tnf\tnf_dist\tnf_polar\tlabels\n')
    for pp_model in list_models:
        p = P()
        if pp_model.endswith('.xyz'):
            p.load_from_xyz(pp_model)
        elif pp_model.endswith('.pma'):
            p.load_from_pma(pp_model)
        else:
            sys.stderr.write(f"This pharmacophore format is not supported. Input file is {pp_model}")

        res = calc_pharm_graph_desc(p)
        if res:
            w.write(os.path.basename(os.path.splitext(pp_model)[0]) + '\t' + '\t'.join(map(str, res)) + '\n')


if __name__ == '__main__':
    entry_point()
