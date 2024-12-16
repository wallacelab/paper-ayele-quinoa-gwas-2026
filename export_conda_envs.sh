#! /bin/bash


# Gather all the conda environments from this project and export their config files for others to recreate
scriptdir=0_Scripts

# First, gather all conda environments from my scripts using their common naming
envs=`grep --no-filename -E 'conda_.+=' $scriptdir/*.sh | sort | uniq | sed -r -e "s|.+=||" -e "s| .+||"`

# Export each environment
for env in $envs; do
    echo "Exporting environment $env"
    conda env export --name $env > $scriptdir/$env.yml
    #break
done

