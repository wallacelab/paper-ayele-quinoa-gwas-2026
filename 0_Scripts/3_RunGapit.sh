#! /bin/bash

# Run GWAS in GAPIT with resampling and permutation analysis

#   Scripts taken from /home/jgwall/Projects/Sorghum/Punurri_aphid_resistance/FarmCPU-RMIP-analysis-v1.3_flowering_with_height_cov/

# Parameters
scriptdir=$1  # Script directory
workdir=$2  # Working directory to output results
genos=$3  # Genotype file from above, modified by coauthor to remove contigs and rename chromosomes
phenos=$4  # Phenotype file
num_pcs=$5   # Number of principal components to include
num_perms=$6  # Number of permutations to run to get null distribution of results
num_resamples=$7 # Number of resamples to run to get RMIPs
resample_fraction=$8 # Resample fraction number
perm_p_cutoff=$9  # P-value cutoff for permutations to count as null "hits"
rmip_cutoff=${10} # RMIP cutoff for a hit to count as real
num_cores=${11} # Number of CPU cores to use


########
# Conda envs
########

# Conda environments
conda_tassel5=tassel-5.2.89

# Make it so can run conda from script
. $(conda info --root)/etc/profile.d/conda.sh 


###########
# Setup
###########

# # Convert genotypes to numeric & generate kinship matrix
# conda activate $conda_tassel5
# run_pipeline.pl -h $genos -NumericalGenotypePlugin -endPlugin -export $workdir/3a_genos.numeric.txt -exportType ReferenceProbability
# run_pipeline.pl -h $genos -KinshipPlugin -endPlugin -export $workdir/3a_kinship.txt
# 
# 
# # Extract genotype coordinates for GAPIT map file
# echo -e "SNP\tChromosome\tPosition" > $workdir/3a_genos.map.txt
# zcat $genos | tail -n +2 | cut -f1,3,4  >>  $workdir/3a_genos.map.txt
# 
# # Combine and check
# Rscript $scriptdir/3a_CombineAndCheckData.r --genofile $workdir/3a_genos.numeric.txt --phenofile $phenos --outprefix $workdir/3a


###########
# Permutations for p-value cutoffs
###########

# Get list of traits
traits=`head -n 1 $workdir/3a_phenos.csv | tr ',' '\n' | tail -n +2 | tr '"' ' '`

# Key files (putting here to make clearer)
phenofile=$workdir/3a_phenos.csv
genofile=$workdir/3a_genos.numeric.rds
mapfile=$workdir/3a_genos.map.txt
kinfile=$workdir/3a_kinship.txt
 
# Loop over traits - TODO - Try with both FarmCPU and with BLINK

for model in FarmCPU BLINK; do

    # Make subdirectories
    permdir=$workdir/3b_Permutations.$model
    resampledir=$workdir/3d_Resampling.$model
    resultdir=$workdir/3e_Results.$model
    if [ ! -e $permdir ] ; then mkdir $permdir; fi
    if [ ! -e $resampledir ] ; then mkdir $resampledir; fi
    if [ ! -e $resultdir ] ; then mkdir $resultdir; fi


    # Loop over traits
    for trait in $traits; do
            
    # # # # #     # Set up Covariates. - Not for this analysis
    # # # # #     case $trait in
    # # # # #         "FLoweringtime") covariates="$cov_height" ;;  
    # # # # #         "Log_2019.2020_PH") covariates="NONE" ;;  
    # # # # #         *) covariates="$cov_flower" ;;  # Just flowering for everything else (=aphid counts and plant height)
    # # # # #     esac
    # # # # #     echo "Running trait $trait with covariates $covariates"
        
    
        # Output directories
        mypermdir=$permdir/$trait
        myresampledir=$resampledir/$trait
        if [ ! -e $mypermdir ]; then mkdir $mypermdir; fi
        if [ ! -e $myresampledir ]; then mkdir $myresampledir; fi
        
        # Run permutations & resamplings
        Rscript $scriptdir/3b_ResampleGapit.r --genofile $genofile --phenofile $phenofile --map $mapfile --kinshipfile $kinfile --trait $trait --num-pcs $num_pcs \
            --permdir $mypermdir --resampledir $myresampledir --num-perms $num_perms --num-resamples $num_resamples --resample-fraction $resample_fraction \
            --model $model --num-cores $num_cores

        # Save actual permutation results with p-values; remove all other intermediate files
        pigz $mypermdir/*/GAPIT*GWAS_Results.$model.$trait.csv
        rm $mypermdir/*/*.pdf $mypermdir/*/*.csv
        
        # Save actual resample results with p-values; remove all other intermediate files
        pigz $myresampledir/*/GAPIT*GWAS_Results.$model.$trait.csv
        rm $myresampledir/*/*.pdf $myresampledir/*/*.csv
         
        
        # Collect permuted p-values
        myperms=$mypermdir/*/GAPIT*GWAS_Results.$model.$trait.csv.gz
        myresamples=$myresampledir/*/GAPIT*GWAS_Results.$model.$trait.csv.gz   
        Rscript $scriptdir/3e_GetRmips.r --perms $myperms --resamples $myresamples --perm-p-cutoff $perm_p_cutoff  --rmip-cutoff $rmip_cutoff \
            -o $resultdir/1e_$trait.p05.rmip05 -t $trait
        
#         break
    done

    #break
done

# TODO: Get linkage among all the RMIP hits and see if any should be combined into single hits?
