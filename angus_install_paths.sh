#!/usr/bin/env bash

## If the user doesn't have a default path to a program, add that to PATH

[[ ":$PATH:" != *":/usr/bin/bowtie2:"* ]] && PATH="/usr/bin/bowtie2:${PATH}"
[[ ":$PATH:" != *":/s/angus/index/common/tools/csa:"* ]] && PATH="/s/angus/index/common/tools/csa:${PATH}"
[[ ":$PATH:" != *":/s/angus/index/common/tools/csa:"* ]] && PATH="/s/angus/index/common/tools/csa:${PATH}"

if [[ -v $genome ]]; echo 'genome="/s/angus/index/databases/bowtie2_indexes/mod_bos_taurus/mod_bos_taurus.fna"' >> ${HOME}/.bashrc; fi
if [[ -v $amr_db ]]; echo 'amr_db="/s/angus/index/databases/MEGARes/megares_database.fasta"' >> ${HOME}/.bashrc; fi
if [[ -v $kraken_db ]]; echo 'kraken_db="/s/angus/index/databases/kraken_databases/Standard_kraken_10.14.db"' >> ${HOME}/.bashrc; fi

