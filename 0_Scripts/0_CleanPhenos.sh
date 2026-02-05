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

# Clean each individually, including outputting QC plots and saving both filtered data and after combining into BLUEs
# Rscript $scriptdir/0b_CleanPhenoData.r -i $workdir/0a_phenos.2024.csv -p $problems2024 --prob-min-row $prob_min_row --prob-max-row $prob_max_row -o $workdir/0b_phenos.2024
Rscript $scriptdir/0b_CleanPhenoData.r -i $workdir/0a_phenos.2025.csv -p $problems2025 --prob-min-row $prob_min_row --prob-max-row $prob_max_row -o $workdir/0b_phenos.2025

# 2025 Tleaf39 looks like a problem, also Tleaf60? E66, VPD66, ETR66, TLeaf66, VPD75, TLeaf75
# TODO - Put 2025 data in alphabetical order of traits and also do problem detection and removal on it
