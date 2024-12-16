#! /bin/bash

# Initialize repo
echo "# proj-quinoa-genotyping-2025" >> README.md
git init
git add README.md
git commit -m "first commit"
git branch -M main
git remote add origin git@github.com:jgwall/proj-quinoa-genotyping-2025.git
git push -u origin main

# Add all script files
find ./ | grep -i -e "\.sh$" -e "\.py$" -e "\.r$" | xargs -d "\n" git add

# Add support files

# Add key intermediate files

# Add output SVGs

# Remove any deleted/untracked files
git add -u

# Commit
timestamp=`date`
git commit -m "$timestamp $1"  # '$1' Includes any extra message typed on command line
git push -u origin main

