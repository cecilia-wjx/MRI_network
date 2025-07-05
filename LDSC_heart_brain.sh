#!/bin/bash

input_file="/cmd_heart_brain.txt"
output_dir="/cmd_heart_brain"
ldsc_script="/ldsc.py"
ref_ld_chr="/eur_w_ld_chr"
w_ld_chr="/eur_w_ld_chr"

# Initialize the row counter
i=1

# Loop through each line in the input file
tail -n +2 "$input_file" | while IFS=' ' read -r file1 file2; do
    # Execute the ldsc.py command
    file2=$(echo "$file2" | sed 's/^[ \t]*//' | tr -d '\r')
    python $ldsc_script \
        --rg "$file1","$file2" \
        --ref-ld-chr "$ref_ld_chr" \
        --w-ld-chr "$w_ld_chr" \
        --out "${output_dir}row${i}"
    
    # Increment the row counter
    i=$((i+1))
done
