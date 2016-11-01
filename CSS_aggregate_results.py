#!/usr/bin/env python

import sys
import glob
import numpy as np

taxa_level = {'D': 0, 'P': 1, 'C': 2, 'O':3, 'F': 4, 'G': 5, 'S': 6}


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
                current_taxon = current_taxon[:taxa_level[entry[3]]]
                current_taxon.append(entry[5])
                print(entry[5], current_taxon)
                taxon_name = '|'.join(current_taxon)
                try:
                    ret[sample_id].setdefault(taxon_name, float(entry[2]))
                except KeyError:
                    ret.setdefault(sample_id, {taxon_name: float(entry[2])})
    return ret


if __name__ == '__main__':
    K = load_kraken_results(sys.argv[1])
    print(K)






