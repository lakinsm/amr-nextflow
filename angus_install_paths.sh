#!/usr/bin/env bash

## If the user doesn't have a default path to a program, add that to PATH

[[ ":$PATH:" != *":/usr/bin/bowtie2:"* ]] && PATH="/usr/bin/bowtie2:${PATH}"
[[ ":$PATH:" != *":/s/angus/index/common/tools/csa:"* ]] && PATH="/s/angus/index/common/tools/csa:${PATH}"

echo "PATH=${PATH}" >> ${HOME}/.bashrc

exit 0

