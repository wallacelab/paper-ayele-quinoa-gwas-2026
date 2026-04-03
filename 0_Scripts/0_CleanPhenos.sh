#! /bin/bash

# Combine phenos and make ready for GWAS

# Arguments
scriptdir=$1 # Folder with accessory scripts in it
workdir=$2  # Working directory for this script
problems2024="$3"
problems2025="$4"
prob_min_row=$5
prob_max_row=$6
priority_traits="$7"  # Traits to prioritize analysis of

##################
# 2024 greenhouse data
##################

# # Standardize phenotype data
# Rscript $scriptdir/0a_StandardizePhenoData.r -g "$workdir/2024_Greenhouse data.csv" -p "$workdir/2025-Physiological traits.xlsx" \
#     -s "$workdir/2025 Seed yield.xlsx" -u "$workdir/2025_3TRT-Plant height _MTD_Yield2.xlsx" -o $workdir/0a_phenos 



# # First pass clean each individually, including outputting QC plots and saving filtered data. Outlier detection is light because of wanting to transform later & re-filter
# # Done outside a loop because each year has specific problem traits need to deal with
# iqr=3
# zscore=5
# Rscript $scriptdir/0b_CleanPhenoData.r -i $workdir/0a_phenos.2024.csv -p $problems2024 --prob-min-row $prob_min_row --prob-max-row $prob_max_row \
#     --iqr-cutoff $iqr --zscore-cutoff $zscore -o $workdir/0b_phenos.2024
# Rscript $scriptdir/0b_CleanPhenoData.r -i $workdir/0a_phenos.2025.csv -p $problems2025 --prob-min-row $prob_min_row --prob-max-row $prob_max_row \
#     --iqr-cutoff $iqr --zscore-cutoff $zscore -o $workdir/0b_phenos.2025
# Rscript $scriptdir/0b_CleanPhenoData.r -i $workdir/0a_phenos.2025_update.csv \
#     --iqr-cutoff $iqr --zscore-cutoff $zscore -o $workdir/0b_phenos.2025_update


# # Boxcox and heavier cleaning. Year-specific problems already done, so can do in a loop
# iqr=1.7
# zscore=3
# for set in 2024 2025 2025_update; do
#     
#     # Boxcox transform
#     Rscript $scriptdir/0c_TransformPhenos.r -i $workdir/0b_phenos.$set.filtered.csv -o $workdir/0c_phenos.$set.filtered
#     
#     # Heavy outlier filtering  
#     Rscript $scriptdir/0b_CleanPhenoData.r -i $workdir/0c_phenos.$set.filtered.transformed.csv \
#         --iqr-cutoff $iqr --zscore-cutoff $zscore -o $workdir/0d_phenos.$set.clean
#         
#     # Format for GAPIT
#     Rscript $scriptdir/0e_FormatPhenosForGapit.r -i $workdir/0d_phenos.$set.clean.filtered.csv -o $workdir/0e_phenos.$set --priority-traits $priority_traits
# done


# Cut out 2024 DTM for separate analysis (easier than trying to hack the pipeline to handle it separately
Rscript -e "x=read.delim('$workdir/0e_phenos.2024.remainder.gapit.txt', check.names=FALSE); tokeep = names(x)[grep(names(x), pattern='^DTM')];" \
    -e "out=x[,c('<Taxa>', tokeep)]; write.table(out, file='$workdir/0e_phenos.2024.DTM.gapit.txt', row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')"

