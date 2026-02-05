#! /bin/bash

# Combine phenos and make ready for GWAS

# Arguments
scriptdir=$1 # Folder with accessory scripts in it
workdir=$2  # Working directory for this script
problems2024="$3"
problems2025="$4"
prob_min_row=$5
prob_max_row=$6

##################
# 2024 greenhouse data
##################

# # Standardize phenotype data
# Rscript $scriptdir/0a_StandardizePhenoData.r -g "$workdir/2024_Greenhouse data.csv" -p "$workdir/2025-Physiological traits.xlsx" -s "$workdir/2025 Seed yield.xlsx" \
#     -o $workdir/0a_phenos 

# TODO - Redo outlier detection to operate on each TRT separately
# TODO - Select target phenos and format for GAPIT?

# # First pass clean each individually, including outputting QC plots and saving filtered data. Outlier detection is light because of wanting to transform later & re-filter
# iqr=3
# zscore=5
# Rscript $scriptdir/0b_CleanPhenoData.r -i $workdir/0a_phenos.2024.csv -p $problems2024 --prob-min-row $prob_min_row --prob-max-row $prob_max_row \
#     --iqr-cutoff $iqr --zscore-cutoff $zscore -o $workdir/0b_phenos.2024
# Rscript $scriptdir/0b_CleanPhenoData.r -i $workdir/0a_phenos.2025.csv -p $problems2025 --prob-min-row $prob_min_row --prob-max-row $prob_max_row \
#     --iqr-cutoff $iqr --zscore-cutoff $zscore -o $workdir/0b_phenos.2025

# # Boxcox and heavier transform. Problems already done, so can do in a loop
# for set in 2024 2025; do
#     Rscript $scriptdir/0c_TransformPhenos.r -i $workdir/0b_phenos.$set.filtered.csv -o $workdir/0c_phenos.$set.filtered
#     
#     # Heavy outlier filtering
#     iqr=1.7
#     zscore=3
#     Rscript $scriptdir/0b_CleanPhenoData.r -i $workdir/0c_phenos.$set.filtered.transformed.csv \
#         --iqr-cutoff $iqr --zscore-cutoff $zscore -o $workdir/0d_phenos.$set.clean
# done




# Boxcox transform everything then re-filter for outliers



# TODO - Remove BLUEs? Nah, just keep them in and don't worry about them.
