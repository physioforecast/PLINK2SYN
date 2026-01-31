#!/bin/bash


syntest_toAv() {
    local gtype="$1"
    local sourcepath="$2"
    local derivativepath="$3"
    local pycodepath="$4"
    local fileroot="1kgp3_50k"
    local pfile="all_hg38"
    
    # Determine the PLINK options based on the gtype
    if [ "$gtype" == "nomiss" ]; then
        local geno_option="--geno 0"
    else
        local geno_option=""
    fi

    mkdir -p "$derivativepath"

    # Run PLINK command with the appropriate options
    plink2 --pfile "${sourcepath}$pfile" vzs \
        --snps-only --max-alleles 2 \
        --make-pfile \
        $geno_option \
        --thin-count 50000 --seed 111 --threads 1 --memory 8000 require \
        --out "${derivativepath}${fileroot}_$gtype"

    # Export to Av format
    plink2 --pfile "${derivativepath}${fileroot}_$gtype" --export Av --out "${derivativepath}${fileroot}_${gtype}_Av"

    # Total number of samples
    total_samples=32017

    # Desired chunk size
    chunk_size=8000

    # Run augment_1kGP3.py with --miss option if gtype is "yesmiss"
    # Loop to create chunks
    start_sample=1
    outputflag=1

    while [ $start_sample -le $total_samples ]; do
        # Calculate end of chunk (make sure it doesn't exceed total_samples)
        end_sample=$((start_sample + chunk_size - 1))
        if [ $end_sample -gt $total_samples ]; then
            end_sample=$total_samples
        fi
        
        # Calculate the actual number of samples in this chunk
        nsamples=$((end_sample - start_sample + 1))

        # Print the current outputflag
        echo "Processing chunk: $outputflag"
    
        
        # Run the Python script for this chunk with appropriate parameters
    # Determine which Python to use
    if command -v python &>/dev/null; then
        PYTHON_CMD=python
    elif command -v python3 &>/dev/null; then
        PYTHON_CMD=python3
    else
        echo "Error: Neither python nor python3 found in PATH." >&2
        exit 1
    fi

    # Run the command
    if [ "$gtype" == "yesmiss" ]; then
        $PYTHON_CMD "${pycodepath}augment_1kGP3.py" \
            --fname "${derivativepath}${fileroot}_${gtype}_Av.traw" \
            --fname2 "${sourcepath}all_hg38.psam" \
            --miss --p 0.1 \
            --nsamples "$nsamples" \
            --outputflag "$outputflag"
    else
        $PYTHON_CMD "${pycodepath}augment_1kGP3.py" \
            --fname "${derivativepath}${fileroot}_${gtype}_Av.traw" \
            --fname2 "${sourcepath}all_hg38.psam" \
            --nsamples "$nsamples" \
            --outputflag "$outputflag"
    fi
        # Update for the next chunk
        start_sample=$((end_sample + 1))
        outputflag=$((outputflag + 1))
    done
}

export -f syntest_toAv