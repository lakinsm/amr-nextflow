#!/usr/bin/env python3

import sys
import glob

amr_level_names = {0: 'Gene', 1: 'Class', 2: 'Mechanism', 3: 'Group'}


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


def amr_aggregate_results(outpath, A, L, N):
    with open(outpath + '/AMR_aggregated_output.csv', 'w') as out:
        out.write('Level,Sample,Name,Count\n')
        for file in glob.glob(outpath + '/amr_output/*.tab'):
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
                    out.write('Gene,{},{},{}\n'.format(
                        sample_id,
                        result[2],
                        result[4]
                    ))
            for name, count in long_results.items():
                out.write('{},{},{},{}\n'.format(
                    N[L[name]],
                    sample_id,
                    name,
                    count
                ))


# def kraken_aggregate_results(outpath)

if __name__ == '__main__':
    A, L = load_megares_annotations(sys.argv[1])
    amr_aggregate_results(sys.argv[2], A, L, amr_level_names)
