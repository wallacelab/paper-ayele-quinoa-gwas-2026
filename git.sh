#! /bin/bash

# # Initialize repo
# echo "# proj-quinoa-genotyping-2025" >> README.md
# git init
# git add README.md
# git commit -m "first commit"
# git branch -M main
# git remote add origin git@github.com:jgwall/proj-quinoa-genotyping-2025.git
# git push -u origin main

# Add all script files
find ./ | grep -i -e "\.sh$" -e "\.py$" -e "\.r$" | xargs -d "\n" git add

# Add support files
git add 0_Scripts/adapters.fa  # Sequencing adapter file for read cleaning
git add 0_Scripts/*.yml # Conda environment exports
git add 0_Scripts/*.csv
git add 0_Phenotypes/0a*.csv # Initial (raw) data
git add 0_Phenotypes/*clean.csv  # Cleaned data after removing outliers and doing boxcox transformation

# Add key intermediate files
git add 3_GWAS/3e_Results*/*.csv
git add 3_GWAS/3f_genome_index.txt

# Add output SVGs


# Remove any deleted/untracked files
git add -u

# Commit
timestamp=`date`
git commit -m "$timestamp $1"  # '$1' Includes any extra message typed on command line
git push -u origin main

