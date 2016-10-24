#!/usr/bin/env bash

## If the user doesn't have a default path to a program, add that to PATH

[[ ":$PATH:" != *":/usr/bin/bowtie2:"* ]] && PATH="/usr/bin/bowtie2:${PATH}"
[[ ":$PATH:" != *":/s/angus/index/common/tools/csa:"* ]] && PATH="/s/angus/index/common/tools/csa:${PATH}"
[[ ":$PATH:" != *":/s/angus/index/common/tools/csa:"* ]] && PATH="/s/angus/index/common/tools/csa:${PATH}"

echo "PATH=${PATH}" >> ${HOME}/.bashrc

if [[ ! $genome ]]; then echo 'genome="/s/angus/index/databases/bowtie2_indexes/mod_bos_taurus/mod_bos_taurus.fna"' >> ${HOME}/.bashrc; fi
if [[ ! $amr_db ]]; then echo 'amr_db="/s/angus/index/databases/MEGARes/megares_database.fasta"' >> ${HOME}/.bashrc; fi
if [[ ! $kraken_db ]]; then echo 'kraken_db="/s/angus/index/databases/kraken_databases/Standard_kraken_10.14.db"' >> ${HOME}/.bashrc; fi

source ~/.bashrc

exit 0

