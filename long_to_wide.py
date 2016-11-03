#!/usr/bin/env python3

import sys


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


def kraken_load_long_css_data(file):
    samples = {}
    labels = {}
    with open(file, 'r') as f:
        data = f.read().split('\n')
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
    for flevel, sdict in S.items():
        local_sample_names = []
        with open(outdir + '/' + flevel + '.csv', 'w') as out:
            for sample in sdict.keys():
                local_sample_names.append(sample)
            out.write(','.join(local_sample_names) + '\n')
            for label in L[flevel]:
                local_counts = []
                out.write(label + ',')
                for sample in local_sample_names:
                    if label in sdict[sample]:
                        local_counts.append(str(sdict[sample][label]))
                    else:
                        local_counts.append(str(0))
                out.write(','.join(local_counts) + '\n')


if __name__ == '__main__':
    S, L = kraken_load_long_data(sys.argv[1])
    output_wide_data(sys.argv[2], S, L)
    del S
    del L
    S, L = kraken_load_long_css_data(sys.argv[3])
    output_wide_data(sys.argv[4], S, L)
    del S
    del L
    S, L = amr_load_long_data(sys.argv[5])
    output_wide_data(sys.argv[6], S, L)

