#!/bin/bash

GWAS_TOOLS=/Xcan/summary-gwas-imputation/src
MRI=/TWAS/brain_data/brain/txt_xcan
DATA=/Xcan/cad_data/data
OUTPUT=/Xcan/mri_output/
METAXCAN=/Xcan/MetaXcan/software

tissues=(
   "Adipose_Visceral_Omentum"
   "Adipose_Subcutaneous"
   "Whole_Blood"
   "Cells_EBV-transformed_lymphocytes"
   "Artery_Coronary"
   "Artery_Aorta" 
   "Artery_Tibial"
   "Brain_Amygdala"
   "Brain_Cerebellum"
   "Brain_Cortex"
   "Brain_Hypothalamus"
   "Brain_Nucleus_accumbens_basal_ganglia"
   "Brain_Anterior_cingulate_cortex_BA24"
   "Brain_Cerebellar_Hemisphere"
   "Brain_Frontal_Cortex_BA9"
   "Brain_Hippocampus"
   "Brain_Putamen_basal_ganglia"
   "Brain_Spinal_cord_cervical_c-1"
   "Brain_Substantia_nigra"
   "Brain_Caudate_basal_ganglia"
   "Heart_Left_Ventricle"
   "Heart_Atrial_Appendage"
   "Liver"
   "Pancreas"
)

count=0
for pheno in "${phenos[@]}"; do
    ((count++))
    
    if [ $((count % 39)) -eq 0 ]; then
        echo "Processing phenotype: ${pheno}"
        
        python "$GWAS_TOOLS/gwas_parsing.py" \
        -gwas_file "$MRI/${pheno}.txt" \
        -liftover "$DATA/liftover/hg19ToHg38.over.chain.gz" \
        -output "$OUTPUT/harmonized_gwas/${pheno}.txt.gz" \
        -snp_reference_metadata "$DATA/reference_panel_1000G/variant_metadata.txt.gz" METADATA \
        -output_column_map rsid variant_id \
        -output_column_map A2 non_effect_allele \
        -output_column_map A1 effect_allele \
        -output_column_map beta effect_size \
        -output_column_map chr chromosome \
        -output_column_map se standard_error \
        -output_column_map p pvalue \
        -output_column_map sample_size sample_size \
        --chromosome_format \
        -output_column_map pos position \
        --insert_value n_cases NA \
        --insert_value frequency NA \
        -output_order variant_id panel_variant_id chromosome position effect_allele non_effect_allele frequency pvalue zscore effect_size standard_error sample_size n_cases

        for chr in {1..22}; do
            for batch in {0..9}; do
                echo "Processing chromosome ${chr}, sub-batch ${batch}"
                output_file="${OUTPUT}/summary_imputation_1000G/${pheno}_chr${chr}_sb${batch}_reg0.1_ff0.01_by_region.txt.gz"
                
                python "$GWAS_TOOLS/gwas_summary_imputation.py" \
                -by_region_file "$DATA/eur_ld.bed.gz" \
                -gwas_file "$OUTPUT/harmonized_gwas/${pheno}.txt.gz" \
                -parquet_genotype "$DATA/reference_panel_1000G/chr${chr}.variants.parquet" \
                -parquet_genotype_metadata "$DATA/reference_panel_1000G/variant_metadata.parquet" \
                -window 100000 \
                -parsimony 7 \
                -chromosome "${chr}" \
                -regularization 0.1 \
                -frequency_filter 0.01 \
                -sub_batches 10 \
                -sub_batch "${batch}" \
                --standardise_dosages \
                -output "${output_file}"
            done
        done

        python "$GWAS_TOOLS/gwas_summary_imputation_postprocess.py" \
        -gwas_file "$OUTPUT/harmonized_gwas/${pheno}.txt.gz" \
        -folder "$OUTPUT/summary_imputation_1000G" \
        -pattern "${pheno}.*" \
        -parsimony 7 \
        -output "$OUTPUT/processed_summary_imputation_1000G/imputed_${pheno}.txt.gz"

        for tissue in "${tissues[@]}"; do
            echo "Processing tissue: ${tissue}"
            
            python "$METAXCAN/SPrediXcan.py" \
            --gwas_file "$OUTPUT/processed_summary_imputation_1000G/imputed_${pheno}.txt.gz" \
            --snp_column panel_variant_id \
            --effect_allele_column effect_allele \
            --non_effect_allele_column non_effect_allele \
            --zscore_column zscore \
            --model_db_path "$DATA/models/eqtl/mashr/mashr_${tissue}.db" \
            --covariance "$DATA/models/eqtl/mashr/mashr_${tissue}.txt.gz" \
            --keep_non_rsid \
            --additional_output \
            --model_db_snp_key varID \
            --throw \
            --output_file "$OUTPUT/spredixcan/eqtl/${pheno}_PM__${tissue}.csv"
        done
    else
        echo "Skipping phenotype ${pheno} (count: ${count})"
    fi
done

