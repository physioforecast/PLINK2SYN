#!/bin/bash
set -e

## Need to create lists of individuals to keep for subsetting.
## To be used for Plink and R matched glm testing sets.
## Need separate sets for nomissing and yesmissing datasets.

# Source config for paths
source ../config.sh 2>/dev/null || source ./config.sh 2>/dev/null || source config.sh

# Validate required variables
if [[ -z "$testdatapath" ]] || [[ -z "$resultsdir" ]]; then
    echo "ERROR: testdatapath and resultsdir must be defined in config.sh"
    echo "Please add:"
    echo "  export testdatapath='test_data/'"
    echo "  export resultsdir='PLINK2SYN_GLM_TESTS_YYYYMMDD/'"
    exit 1
fi

datapath="${testdatapath}${resultsdir}"

fnames_=(
  "1kgp3_50k_nomiss_Av_nonintdose_recode_varIDs.psam"
  "1kgp3_50k_yesmiss_Av_nonintdose_recode_varIDs.psam"
) 

## We can grab n rows, column 1 (#IID) from psam.
thin=(1000 32017 32000)

for fname_only in "${fnames_[@]}"; do   # <-- fixed: fnames_ instead of fnames
  infile="${datapath}${fname_only}"

  if [[ ! -f "$infile" ]]; then
    echo "⚠️  Skipping: $infile not found."
    continue
  fi

  total_lines=$(wc -l < "$infile")
  if (( total_lines <= 1 )); then
    echo "⚠️  Skipping empty or header-only file: $infile"
    continue
  fi

  data_rows=$(( total_lines - 1 ))

  for n in "${thin[@]}"; do
    if (( n > data_rows )); then
      echo "⚠️  Requested $n rows but $infile has only $data_rows — using $data_rows instead."
      used_n=$data_rows
    else
      used_n=$n
    fi

    outname="${datapath}$(basename "${infile%.psam}")_subset_${used_n}.keep"
    echo "Creating subset: $outname (header + $used_n data rows) from $infile"

    # Header (1 row)
    head -n 1 "$infile" | cut -f1 > "$outname"

    # First column of next used_n rows
    tail -n +2 "$infile" | head -n "$used_n" | cut -f1 >> "$outname"

    # Verify line count
    lines=$(wc -l < "$outname")
    if [[ "$lines" -eq $((used_n + 1)) ]]; then
      echo "✅ Verified: $outname has $lines lines (1 header + $used_n data rows)"
    else
      echo "❌ Mismatch: $outname has $lines lines, expected $((used_n + 1))"
    fi
  done
done

