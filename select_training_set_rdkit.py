#!/usr/bin/env python3
# author          : Alina Kutlushina
# date            : 01.05.2018
# license         : BSD-3
#==============================================================================

import os
import sys
import argparse
from collections import defaultdict
from rdkit import Chem, DataStructs
from rdkit.ML.Cluster import Butina
from rdkit.Chem import AllChem, ChemicalFeatures
from rdkit.Chem.Pharm2D import Generate
from rdkit.Chem.Pharm2D.SigFactory import SigFactory


def read_file(fname, fcfp4, fdef_fname):
    if not fcfp4:
        featFactory = ChemicalFeatures.BuildFeatureFactory(fdef_fname)
        sigFactory = SigFactory(featFactory, minPointCount=2, maxPointCount=3, trianglePruneBins=False)
        sigFactory.SetBins([(0, 2), (2, 5), (5, 8)])
        sigFactory.Init()
    d = defaultdict(list)
    with open(fname) as f:
        for row in f:
            smiles, ids, aff = row.strip().split('\t')
            if smiles is not None:
                mol = Chem.MolFromSmiles(smiles)
                d['mol_name'].append(ids)
                d['smiles'].append(smiles)
                if fcfp4:
                    d['fingerprint'].append(AllChem.GetMorganFingerprint(mol, 2, useFeatures=True))
                else:
                    d['fingerprint'].append(Generate.Gen2DFingerprint(mol, sigFactory))
    return d


def gen_cluster_subset_algButina(fps, cutoff):
    dists = []
    for i, fp in enumerate(fps):
        distance_matrix = DataStructs.BulkTanimotoSimilarity(fps[i], fps[:i])
        dists.extend([1 - x for x in distance_matrix])
    cs = Butina.ClusterData(dists, len(fps), cutoff, isDistData=True)
    return cs  # returns tuple of tuples with sequential numbers of compounds in each cluster


def save_cluster_stat(cs_index, len_act, clust_stat):
    for i, cluster in enumerate(cs_index):
        i_act = 0
        for el in cluster:
            if el in range(0, len_act+1):
               i_act += 1
        #print('cluster №%i, cluster length %i, share of active %.2f' % (i, len(cluster), i_act/len(cluster)))
        #print(cluster, '\n')
        clust_stat.write('cluster №%i, cluster length %i, share of active %.2f \n' % (i, len(cluster), i_act/len(cluster)))


def diff_binding_mode(cs, mol_names, smiles, len_act, inact_centroids, min_num):
    # mol_names contains actives and then inactives
    # therefore len_act is equal to the number of actives to select them from the list of mol names
    inact_centroids = set(inact_centroids)  # set of tuple of tuples
    ts_full = []
    for i, c in enumerate(cs):
        if len(c) > min_num:
            ts_mol_name_act = []
            ts_mol_name_inact = []
            for x in c:
                if x in range(0, len_act + 1) and len(ts_mol_name_act) < min_num:
                    ts_mol_name_act.append((mol_names[x], smiles[x]))  # list of tuples
                elif x not in range(0, len_act + 1) and len(ts_mol_name_inact) < min_num:
                    ts_mol_name_inact.append((mol_names[x], smiles[x]))
            ts_mol_name_inact = set(ts_mol_name_inact) | inact_centroids
            if len(ts_mol_name_act) == min_num:
                ts_full.append((i, tuple(sorted(ts_mol_name_act)), tuple(sorted(ts_mol_name_inact))))  # list of two tuples of tuples
    return tuple(ts_full)


def get_centroids(cs, d_msf, num):
    return tuple(sorted((d_msf['mol_name'][x[0]], d_msf['smiles'][x[0]]) for x in cs if len(x) >= num))


def main(in_fname_act, in_fname_inact, output,
         fdef_fname, make_clust, fcfp4, clust_stat, threshold_clust, clust_size, max_nact_trainset):

    d_msf_act = read_file(in_fname_act, fcfp4, fdef_fname)  # defaltdict(list) keys: mol_name, SMILES, fingerprint
    d_msf_inact = read_file(in_fname_inact, fcfp4, fdef_fname)  # defaltdict(list) keys: mol_name, SMILES, fingerprint

    ftrainset = []

    if make_clust:
        cs = gen_cluster_subset_algButina(d_msf_act['fingerprint'] + d_msf_inact['fingerprint'], threshold_clust)
        cs_inact = gen_cluster_subset_algButina(d_msf_inact['fingerprint'], threshold_clust)
        inact_centroids = get_centroids(cs_inact, d_msf_inact, clust_size)  # tuple of tuples with mol names and their SMILES
        if clust_stat:
            save_cluster_stat(cs, len(d_msf_act['mol_name']), clust_stat)
        ts_full = diff_binding_mode(cs, d_msf_act['mol_name'] + d_msf_inact['mol_name'],
                                    d_msf_act['smiles'] + d_msf_inact['smiles'],
                                    len(d_msf_act['mol_name']), inact_centroids, max_nact_trainset)  # tuple of lists of two tuples of tuples
        for i, act_ts, inact_ts in ts_full:
            ftrainset.append([os.path.join(output, 'active_tr%i.csv' % (i)),
                              os.path.join(output, 'inactive_tr%i.csv' % (i))])
            with open(os.path.join(output, 'active_tr%i.csv' % (i)), 'wt') as f:
                f.write('\n'.join('{}\t{}'.format(smiles, mol_name) for mol_name, smiles in act_ts))
            with open(os.path.join(output, 'inactive_tr%i.csv' % (i)), 'wt') as f:
                f.write('\n'.join('{}\t{}'.format(smiles, mol_name) for mol_name, smiles in inact_ts))

    else:
        ftrainset.append([os.path.join(output, 'active_centroid.csv'),
                          os.path.join(output, 'inactive_centroid.csv')])
        # process actives
        cs = gen_cluster_subset_algButina(d_msf_act['fingerprint'], threshold_clust)
        centroids = get_centroids(cs, d_msf_act, clust_size)
        with open(os.path.join(output, 'active_centroid.csv'), 'wt') as f:
            f.write('\n'.join('{}\t{}'.format(smiles, mol_name) for mol_name, smiles in centroids))

        # process inactives
        cs = gen_cluster_subset_algButina(d_msf_inact['fingerprint'], threshold_clust)
        centroids = get_centroids(cs, d_msf_inact, clust_size)
        with open(os.path.join(output, 'inactive_centroid.csv'), 'wt') as f:
            f.write('\n'.join('{}\t{}'.format(smiles, mol_name) for mol_name, smiles in centroids))
    return ftrainset


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='select compounds for training set',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-a', '--in_act', metavar='input_active.smi', required=True,
                        help='input SMILES file name with active compounds.')
    parser.add_argument('-i', '--in_inact', metavar='input_inactive.smi', required=True,
                        help='input SMILES file name with inactive compounds.')
    parser.add_argument('-o', '--output', metavar='output/path', default=None,
                        help='output path. The folder where will be saved a training set.')
    parser.add_argument('-clust', '--make_clust', action='store_true', default=False,
                        help='if set training sets will be created for separate clusters, '
                             'otherwise only one training set will be created.')
    parser.add_argument('-f', '--rdkit_fdef', metavar='smarts.fdef', required=False, default=None,
                        help='fdef-file with pharmacophore feature definition.')
    parser.add_argument('--fcfp4', action='store_true', default=False,
                        help='if set FCFP4 fingerprints will be used for compound selection, '
                             'otherwise pharmacophore fingerprints will be used based on feature '
                             'definitions provided by --rdkit_fdef argument.')
    parser.add_argument('-s', '--cluster_stat', default=None,
                        help='if designate path to file than save cluster statistics')
    parser.add_argument('-t', '--threshold_clust', default=0.4,
                        help='treshold for сlustering data by Butina algorithm')
    parser.add_argument('-clz', '--clust_size', default=5,
                        help='minimum cluster size from extract centroinds for training set')
    parser.add_argument('-m', '--max_act_ts', default=5,
                        help='maximum number of active compounds for training set')

    args = vars(parser.parse_args())
    for o, v in args.items():
        if o == "in_act": in_fname_act = v
        if o == "in_inact": in_fname_inact = v
        if o == "output": output = v
        if o == "rdkit_fdef": fdef_fname = v
        if o == "make_clust": make_clust = v
        if o == "fcfp4": fcfp4 = v
        if o == "cluster_stat": clust_stat = v
        if o == "threshold_clust": threshold_clust = float(v)
        if o == "clust_size": clust_size = int(v)
        if o == "max_act_ts": max_nact_trainset = int(v)

    if not fcfp4 and not fdef_fname:
        sys.stderr.write('At least one argument fcfp4 or rdkit_fdef should be set to a non default value.')
        exit()

    if not output:
        output = os.path.join(os.path.split(os.path.dirname(os.path.abspath(in_fname_act)))[0], 'trainset')
        if not os.path.isdir(output):
            os.makedirs(output)

    main(in_fname_act=in_fname_act,
         in_fname_inact=in_fname_inact,
         output=output,
         fdef_fname=fdef_fname,
         make_clust=make_clust,
         fcfp4=fcfp4,
         clust_stat=clust_stat,
         threshold_clust=threshold_clust,
         clust_size=clust_size,
         max_nact_trainset=max_nact_trainset)
