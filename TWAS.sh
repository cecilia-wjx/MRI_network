#!/bin/bash

sumstats_dir="/dir_path"
weights_dir="/WEIGHT_ALL"
out_dir="/output_dir"
ld_ref="/1000G.EUR."

# Iterate through sumstats files
for sumstats_file in "$sumstats_dir"/*.sumstats; do
    for weights_file in "$weights_dir"/*.pos; do
        for chr in {1..22}; do
            # Extract base names for output clarity
            sumstats_base=$(basename "$sumstats_file" .sumstats)
            weights_base=$(basename "$weights_file" .pos)

            # Define output file name
            output_file="$out_dir/${sumstats_base}_${weights_base}_chr${chr}.dat"

            # Run the R script for FUSION.assoc_test
            Rscript /FUSION.assoc_test.R \
            --sumstats "$sumstats_file" \
            --weights "$weights_file" \
            --weights_dir "$weights_dir" \
            --ref_ld_chr "${ld_ref}" \
            --chr "$chr" \
            --out "$output_file"
        done
    done
done
