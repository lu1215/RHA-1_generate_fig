#!/bin/bash
# -------------------

BASE_DIR=$(dirname $0)/..
task=${BASE_DIR}/../src/srna_analyze.py
config_path=${BASE_DIR}/config

python $task --config $config_path --density --metagene --boundary --codon --fold --scatter
