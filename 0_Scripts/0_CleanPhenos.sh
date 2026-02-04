#! /bin/bash

# Combine phenos and make ready for GWAS

# Arguments
scriptdir=$1 # Folder with accessory scripts in it
workdir=$2  # Working directory for this script
problems2024="$3"

##################
# 2024 greenhouse data
##################

# TODO - Redo these scripts into a "standardize" and "diagnose" script since are 80% identical

Rscript $scriptdir/0a_Clean2024Data.r -i "$workdir/2024_Greenhouse data.csv" -o $workdir/0a_2024_phenos -p $problems2024
# Rscript $scriptdir/0a_Clean2025Data.r -i "$workdir/2025-Physiological traits.xlsx" -s "$phenodir/2025 Seed yield.xlsx" -o $phenodir/0a_2025_phenos

# 2025 Tleaf39 looks like a problem, also Tleaf60? E66, VPD66, ETR66, TLeaf66, VPD75, TLeaf75
# TODO - Put 2025 data in alphabetical order of traits and also do problem detection and removal on it
