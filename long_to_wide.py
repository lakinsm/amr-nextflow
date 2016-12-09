#!/usr/bin/env python3

## sys.argv:
#   1. kraken_raw_output.csv
#   2. nextflow_output folder
#   3. kraken_count_matrices folder
#   4. AMR_raw_output.csv
#   5. amr_count_matrices folder

import sys
import glob
import numpy as np

taxa_level = {'D': 0, 'P': 1, 'C': 2, 'O': 3, 'F': 4, 'G': 5, 'S': 6}
taxa_level_names = {1: 'Domain', 2: 'Phylum', 3: 'Class', 4: 'Order',
                    5: 'Family', 6: 'Genus', 7: 'Species', 8: 'Unclassified'}
amr_level_names = {0: 'Class', 1: 'Mechanism', 2: 'Group'}


def kraken_load_long_data(file):
    samples = {}
    labels = {}
    with open(file, 'r') as f:
        data = f.read().split('\n')
        for entry in data:
            if not entry:
                continue
            entry = entry.split(',')
            try:
                samples[entry[1]][entry[0]][entry[6]] = float(entry[2])
            except KeyError:
                try:
                    samples[entry[1]][entry[0]].setdefault(entry[6], float(entry[2]))
                except KeyError:
                    try:
                        samples[entry[1]].setdefault(entry[0], {entry[6]: float(entry[2])})
                    except KeyError:
                        samples.setdefault(entry[1], {entry[0]: {entry[6]: float(entry[2])}})
            try:
                if entry[6] not in labels[entry[1]]:
                    labels[entry[1]] += (entry[6],)
            except KeyError:
                labels.setdefault(entry[1], (entry[6],))
    return samples, labels


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


def kraken_load_analytic_data(dirpath):
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
    print(ret)
    return dict_to_matrix(ret)


def amr_load_long_data(file):
    samples = {}
    labels = {}
    with open(file, 'r') as f:
        data = f.read().split('\n')[1:]
        for entry in data:
            if not entry:
                continue
            entry = entry.split(',')
            try:
                samples[entry[0]][entry[1]][entry[2]] = float(entry[3])
            except KeyError:
                try:
                    samples[entry[0]][entry[1]].setdefault(entry[2], float(entry[3]))
                except KeyError:
                    try:
                        samples[entry[0]].setdefault(entry[1], {entry[2]: float(entry[3])})
                    except KeyError:
                        samples.setdefault(entry[0], {entry[1]: {entry[2]: float(entry[3])}})
            try:
                if entry[2] not in labels[entry[0]]:
                    labels[entry[0]] += (entry[2],)
            except KeyError:
                labels.setdefault(entry[0], (entry[2],))
    return samples, labels


def output_wide_data(outdir, S, L):
    with open(outdir + '/../AMR_analytic_matrix.csv', 'w') as amr:
        for flevel, sdict in S.items():
            local_sample_names = []
            with open(outdir + '/' + flevel + '.csv', 'w') as out:
                for sample in sdict.keys():
                    local_sample_names.append(sample)
                out.write(','.join(local_sample_names) + '\n')
                if flevel == 'Gene':
                    amr.write(','.join(local_sample_names) + '\n')
                for label in L[flevel]:
                    local_counts = []
                    if flevel == 'Gene':
                        amr.write(label + ',')
                    out.write(label + ',')
                    for sample in local_sample_names:
                        if label in sdict[sample]:
                            local_counts.append(str(sdict[sample][label]))
                        else:
                            local_counts.append(str(0))
                    if flevel == 'Gene':
                        amr.write(','.join(local_counts) + '\n')
                    out.write(','.join(local_counts) + '\n')


def output_kraken_analytic_data(outdir, M, m_names, n_names):
    with open(outdir + '/kraken_analytic_matrix.csv', 'w') as out:
        out.write(','.join(n_names) + '\n')
        for i, row in enumerate(M):
            out.write('{},'.format(m_names[i].replace(',', '_')))
            out.write(','.join([str(x) for x in row]) + '\n')


if __name__ == '__main__':
    S, L = kraken_load_long_data(sys.argv[1])
    K, m, n = kraken_load_analytic_data(sys.argv[2])
    output_kraken_analytic_data(sys.argv[2], K, m, n)
    output_wide_data(sys.argv[3], S, L)
    print(K, m, n)
    del S
    del L
    S, L = amr_load_long_data(sys.argv[4])
    output_wide_data(sys.argv[5], S, L)

