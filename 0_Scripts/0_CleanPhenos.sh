#! /bin/bash

# Combine phenos and make ready for GWAS

# Arguments
scriptdir=$1 # Folder with accessory scripts in it
workdir=$2  # Working directory for this script
problems2024="$3"

##################
# 2024 greenhouse data
##################

Rscript $scriptdir/0a_Clean2024Data.r -i "$workdir/2024_Greenhouse data.csv" -o $workdir/0a_2024_phenos -p $problems2024
# Rscript $scriptdir/0a_Clean2025Data.r -i "$workdir/2025-Physiological traits.xlsx" -s "$phenodir/2025 Seed yield.xlsx" -o $phenodir/0a_2025_phenos
