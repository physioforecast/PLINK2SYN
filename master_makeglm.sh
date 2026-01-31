#!/bin/bash

source config.sh

pycodepath="./code/"
sourcepath="./sourcedata/"

## Make GLM - generate phenotypes with known genetic effects
source PLINK2SYN-venv/bin/activate

# Generate phenotypes for both scenarios (uncomment as needed)
python ${pycodepath}/gen_gwas.py "${derivativepath}1kgp3_50k_nomiss_Av_nonintdose_combined.traw"
python ${pycodepath}/gen_gwas.py "${derivativepath}1kgp3_50k_yesmiss_Av_nonintdose_combined.traw"

