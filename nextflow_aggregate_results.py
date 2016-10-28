#!/usr/bin/env python3

import sys

amr_level_names = {0: 'Gene', 1: 'Class', 2: 'Mechanism', 3: 'Group'}
taxa_level_names = {'D': 'Domain', 'P': 'Phylum', 'C': 'Class', 'O': 'Order',
                    'F': 'Family', 'G': 'Genus', 'S': 'Species', 'U': 'Unclassified'}


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


def amr_aggregate_results(file, A, L, N):
    long_results = {}
    sample_id = file.split('/')[-1].replace('_coverage_sampler_amr.tab', '')
    with open(file, 'r') as f:
        data = f.read().split('\n')[1:]
        for result in data:
            if not result:
                continue
            result = result.split('\t')
            for l in A[result[2]]:
                if l not in long_results:
                    long_results[l] = float(result[4])
                else:
                    long_results[l] += float(result[4])
    for name, count in long_results.items():
        sys.stdout.write('{},{},{},{}\n'.format(
            N[L[name]],
            sample_id,
            name,
            count
        ))


def kraken_aggregate_results(file, T):
    sample_id = file.split('/')[-1].replace('_kraken.report', '')
    with open(file, 'r') as f:
        data = f.read().split('\n')
        for result in data:
            if not result:
                continue
            result = result.split()
            if result[3] == '-':
                continue
            sys.stdout.write('{},{},{},{},{},{},{}\n'.format(
                sample_id,
                T[result[3]],
                result[1],
                result[2],
                result[0],
                result[4],
                ' '.join(result[5:])
            ))


if __name__ == '__main__':
    if sys.argv[1] == 'AMR':
        A, L = load_megares_annotations(sys.argv[2])
        amr_aggregate_results(sys.argv[3], A, L, amr_level_names)
    elif sys.argv[1] == "kraken":
        kraken_aggregate_results(sys.argv[2], taxa_level_names)
    else:
        sys.exit("First argument was neither AMR nor kraken\n")
