#!/usr/bin/env python

import sys
import glob
import numpy as np

np.set_printoptions(suppress=True)
taxa_level = {'D': 0, 'P': 1, 'C': 2, 'O': 3, 'F': 4, 'G': 5, 'S': 6}


def dict_to_matrix(D):
    ncol = len(D.keys())
    unique_nodes = []
    for sample, tdict in D.items():
        for taxon in tdict.keys():
            if taxon not in unique_nodes:
                unique_nodes.append(taxon)
    nrow = len(unique_nodes)
    ret = np.zeros((nrow, ncol), dtype=np.float)
    print(nrow, ncol, len(unique_nodes))
    for j, (sample, tdict) in enumerate(D.items()):
        for i, taxon in enumerate(unique_nodes):
            if taxon in tdict:
                ret[i, j] = np.float(tdict[taxon])
    return ret


def load_kraken_results(dirpath):
    ret = {}
    for file in glob.glob(dirpath + '/kraken_output/*.report'):
        sample_id = file.split('/')[-1].replace('_kraken.report', '')
        with open(file, 'r') as f:
            data = f.read().split('\n')
            current_taxon = []
            for entry in data:
                if not entry:
                    continue
                entry = entry.split()
                if entry[3] in ('-', 'U'):
                    continue
                node_name = ' '.join(entry[5:])
                current_taxon = current_taxon[:taxa_level[entry[3]]]
                current_taxon.append(node_name)
                if float(entry[2]) == 0:
                    continue
                taxon_name = '|'.join(current_taxon)
                try:
                    ret[sample_id].setdefault(taxon_name, float(entry[2]))
                except KeyError:
                    ret.setdefault(sample_id, {taxon_name: float(entry[2])})
    return ret


def calculate_css_norm_factors(K):
    M = dict_to_matrix(K)
    sample_sums = np.array(M.sum(axis=0))
    print(sample_sums)
    qlj = np.zeros(M.shape[1])  # A vector of quantile values for each sample
    l_current = float(0.01)  # The starting point for choice of lth quantile
    while l_current < 1:  # This will change to the stopping condition later
        for j, sample in enumerate(M.T):
            positives = np.sort(np.array(sample[sample > 0]))
            l_iter = 0
            s_lj = positives[l_iter]
            while (float(s_lj) / sample_sums[j]) < l_current:
                l_iter += 1
                s_lj = positives[l_iter]
            qlj[j] = s_lj
        l_current += 0.01
    print(qlj)


if __name__ == '__main__':
    K = load_kraken_results(sys.argv[1])
    calculate_css_norm_factors(K)






