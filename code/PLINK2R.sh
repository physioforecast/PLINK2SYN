#!/bin/bash

set -e

# Source config for paths
source ../config.sh 2>/dev/null || source ./config.sh 2>/dev/null || source config.sh

# Validate required variables
if [[ -z "$testdatapath" ]]; then
    echo "ERROR: testdatapath must be defined in config.sh"
    echo "Please add: export testdatapath='test_data/'"
    exit 1
fi

fname="${testdatapath}1kgp3_50k_nomiss_Av_nonintdose"

plink2 \
    --pfile "$fname" \
    --make-pgen \
    --set-all-var-ids @:#\$r,\$a \
    --out "${fname}_recode_varIDs"

plink2 \
  --pfile "${fname}_recode_varIDs" \
  --export A \
  --out "${fname}_recode_varIDs_A"

fname="${testdatapath}1kgp3_50k_yesmiss_Av_nonintdose"

plink2 \
    --pfile "$fname" \
    --make-pgen \
    --set-all-var-ids @:#\$r,\$a \
    --out "${fname}_recode_varIDs"

plink2 \
  --pfile "${fname}_recode_varIDs" \
  --export A \
  --out "${fname}_recode_varIDs_A"