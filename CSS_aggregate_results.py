#!/usr/bin/env python

import sys
import glob
import numpy as np

np.set_printoptions(suppress=True)
taxa_level = {'D': 0, 'P': 1, 'C': 2, 'O': 3, 'F': 4, 'G': 5, 'S': 6}
taxa_level_names = {1: 'Domain', 2: 'Phylum', 3: 'Class', 4: 'Order',
                    5: 'Family', 6: 'Genus', 7: 'Species', 8: 'Unclassified'}
amr_level_names = {0: 'Class', 1: 'Mechanism', 2: 'Group'}


def dict_to_matrix(D):
    ncol = len(D.keys())
    unique_nodes = []
    samples = []
    for sample, tdict in D.items():
        for taxon in tdict.keys():
            if taxon not in unique_nodes:
                unique_nodes.append(taxon)
    nrow = len(unique_nodes)
    ret = np.zeros((nrow, ncol), dtype=np.float)
    for j, (sample, tdict) in enumerate(D.items()):
        samples.append(sample)
        for i, taxon in enumerate(unique_nodes):
            if taxon in tdict:
                ret[i, j] = np.float(tdict[taxon])
    return ret, unique_nodes, samples


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


def load_megares_annotations(file):
    annot = {}
    levels = {}
    with open(file, 'r') as f:
        data = f.read().split('\n')[1:]
        for entry in data:
            entry = entry.split(',')
            annot.setdefault(entry[0], entry[1:])
            for e, level in enumerate(entry):
                levels.setdefault(level, e)
    return annot, levels


def load_amr_results(dirpath):
    ret = {}
    for file in glob.glob(dirpath + 'amr_output/*.tab'):
        sample_id = file.split('/')[-1].replace('_coverage_sampler_amr.tab', '')
        with open(file, 'r') as f:
            data = f.read().split('\n')[1:]
            for entry in data:
                if not entry:
                    continue
                entry = entry.split('\t')
                try:
                    ret[sample_id].setdefault(entry[2], float(entry[4]))
                except KeyError:
                    ret.setdefault(sample_id, {entry[2]: float(entry[4])})
    return ret


def calculate_css_norm_factors(D):
    M, m_names, n_names = dict_to_matrix(D)
    sample_sums = np.array(M.sum(axis=0))
    qlj = np.zeros(M.shape[1])  # A vector of quantile values for each sample
    last_qlj_idx = [0] * M.shape[1]  # Keep track of where we left off last time
    l_current = float(0.01)  # The starting point for choice of lth quantile
    d_l_previous = np.zeros(M.shape[1])  # Initialize array for calculating stopping condition
    l_chosen = 0
    while l_current < 1:
        # Calculate slj for each sample j
        for j, sample in enumerate(M.T):
            positives = np.sort(np.array(sample[sample > 0]))
            l_iter = last_qlj_idx[j]
            try:
                slj = positives[l_iter]
                while (float(slj) / sample_sums[j]) < l_current:
                    l_iter += 1
                    try:
                        slj = positives[l_iter]
                    except IndexError:
                        slj = positives[-1]
                        break
            except IndexError:
                slj = positives[-1]
            last_qlj_idx[j] = l_iter
            qlj[j] = slj
        # Calculate stopping criterion
        qlmed = np.median(qlj)
        d_l_current = np.median(np.abs(qlj - qlmed))
        if l_current == 0.01:
            d_l_previous = d_l_current
        else:
            if (d_l_current - d_l_previous) >= d_l_previous:
                l_chosen = l_current - 0.01
                if l_chosen < 0.5:
                    l_chosen = 0.5  # Use 0.5 as lowest value, recommended in literature
                break
        # Try next choice of l
        l_current += 0.01
    qlj_ret = np.zeros(M.shape[1])
    for j, sample in enumerate(M.T):
        positives = np.sort(np.array(sample[sample > 0]))
        l_iter = 0
        slj = positives[l_iter]
        while (float(slj) / sample_sums[j]) < l_chosen:
            try:
                slj = positives[l_iter]
            except IndexError:
                slj = positives[-1]
                break
            l_iter += 1
        qlj_ret[j] = slj
    return M, l_chosen, qlj_ret, m_names, n_names


def normalize_matrix(M, N, slj):
    return (M / slj) * N


def kraken_aggregate_and_output(dirpath, M, m_names, n_names):
    ret = {}
    for i, row in enumerate(M):
        feature_name = m_names[i]
        while feature_name:
            try:
                ret[feature_name] += M[i]
            except KeyError:
                ret.setdefault(feature_name, M[i])
            feature_name = feature_name[:-1]
    with open(dirpath + '/kraken_css_normalized_output.csv', 'w') as out:
        for feature, row in ret.items():
            for j, count in enumerate(row):
                out.write('{},{},{},{}\n'.format(
                    taxa_level_names[len(feature.split('|'))],
                    n_names[j],
                    feature,
                    count
                ))


def amr_aggregate_and_output(dirpath, M, m_names, n_names, A, L):
    ret = {}
    level_dict = {}
    for i, row in enumerate(M):
        feature_name = m_names[i]
        ret.setdefault(feature_name, M[i])
        level_dict.setdefault(feature_name, 'Gene')
        for e, term in enumerate(A[feature_name]):
            level_dict.setdefault(term, amr_level_names[e])
            try:
                ret[term] += M[i]
            except KeyError:
                ret.setdefault(term, M[i])
    with open(dirpath + '/amr_css_normalized_output.csv', 'w') as out:
        for feature, row in ret.items():
            for j, count in enumerate(row):
                out.write('{},{},{},{}\n'.format(
                    level_dict[feature],
                    n_names[j],
                    feature,
                    count
                ))

if __name__ == '__main__':
    # Kraken data
    K = load_kraken_results(sys.argv[1])
    K_m, N, slj, m_names, s_names = calculate_css_norm_factors(K)
    M_norm = normalize_matrix(K_m, N, slj)
    kraken_aggregate_and_output(sys.argv[1], M_norm, m_names, s_names)
    del M_norm, m_names, s_names

    # AMR data
    AMR = load_amr_results(sys.argv[1])
    A, L = load_megares_annotations(sys.argv[2])
    AMR_m, N, slj, m_names, s_names = calculate_css_norm_factors(AMR)
    M_norm = normalize_matrix(AMR_m, N, slj)
    amr_aggregate_and_output(sys.argv[1], M_norm, m_names, s_names, A, L)






