#!/bin/bash

# Source config for derivativepath
source ../config.sh 2>/dev/null || source ./config.sh 2>/dev/null || true

sourcepath="./"
derivativepath="${derivativepath:-./tmp/}"  # Use config.sh or default to ./tmp/
pycodepath="./"
gtype="nomiss"
fileroot="1kgp3_50k"
pfile="all_hg38"
    
# Array of file names
file_list=("$@")  # All arguments passed to the script are treated as files

# Create an empty temporary file to store combined content
output_file=${derivativepath}${fileroot}_${gtype}_Av_nonintdose_combined.psam

# Loop through each file
for i in "${!file_list[@]}"; do
  file="${file_list[$i]}"
  
  # Check if it is the first file
  if [ "$i" -eq 0 ]; then
    # For the first file, append all rows (including header)
    cat "${derivativepath}${file}" > "$output_file"
  else
    # For subsequent files, skip the first row (header) and append
    tail -n +2 "${derivativepath}${file}" >> "$output_file"
  fi
done


plink2 --import-dosage ${derivativepath}1kgp3_50k_${gtype}_Av_nonintdose_combined.traw skip0=1 skip1=2 id-delim=_ chr-col-num=1 pos-col-num=4 ref-first --out ${derivativepath}1kgp3_50k_${gtype}_Av_nonintdose --psam ${derivativepath}1kgp3_50k_${gtype}_Av_nonintdose_combined.psam

# plink2 --import-dosage ${derivativepath}1kgp3_50k_${gtype}_Av_nonintdose_5.traw skip0=1 skip1=2 id-delim=_ chr-col-num=1 pos-col-num=4 ref-first --out ${derivativepath}1kgp3_50k_${gtype}_Av_nonintdose --psam ${derivativepath}1kgp3_50k_${gtype}_Av_nonintdose_5.psam
