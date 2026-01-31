#!/bin/bash
# Uses PLINK v2.0.0-a.6.1 M1 (14 Nov 2024) and Python - need to be in path
# Ensure script stops if there is an error and handles pipes correctly
set -exo pipefail

# Source config for derivativepath
source ../config.sh 2>/dev/null || source ./config.sh 2>/dev/null || true

# Function to process PLINK commands
syntest_toAv_combinechunks() {
    local sourcepath="./"
    local derivativepath="${derivativepath:-./tmp/}"  # Use config.sh or default to ./tmp/
    local pycodepath="./"
    local gtype="$1"
    local fileroot="1kgp3_50k"
    local pfile="all_hg38"
    
    # Total number of samples
    total_samples=32017

    # Desired chunk size
    chunk_size=8000

    # Calculate number of chunks (rounding up)
    chunks=$(( (total_samples + chunk_size - 1) / chunk_size ))

    echo "Total chunks: $chunks"

    # Loop over the number of chunks
    for chunk in $(seq 2 $chunks); do
        cut -f7- ${derivativepath}${fileroot}_${gtype}_Av_nonintdose_${chunk}.traw > ${derivativepath}${fileroot}_${gtype}_Av_nonintdose_${chunk}_cut.traw
    done
    
    # Combine all chunk files into one
    paste ${derivativepath}${fileroot}_${gtype}_Av_nonintdose_1.traw $(for chunk in $(seq 2 $chunks); do echo -n "${derivativepath}${fileroot}_${gtype}_Av_nonintdose_${chunk}_cut.traw "; done) > ${derivativepath}${fileroot}_${gtype}_Av_nonintdose_combined.traw


}

# Export function for sourcing
export -f syntest_toAv_combinechunks

# Only run if executed directly (not sourced)
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    syntest_toAv_combinechunks "${1:-nomiss}"
fi
