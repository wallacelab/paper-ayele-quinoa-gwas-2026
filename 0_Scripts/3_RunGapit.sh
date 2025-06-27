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


# Make subdirectories
permdir=$workdir/3b_Permutations
resampledir=$workdir/3d_Resampling
resultdir=$workdir/3e_Results
if [ ! -e $permdir ] ; then mkdir $permdir; fi
if [ ! -e $resampledir ] ; then mkdir $resampledir; fi
if [ ! -e $resultdir ] ; then mkdir $resultdir; fi


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

# # Convert genotypes to numeric
# conda activate $conda_tassel5
# run_pipeline.pl -h $genos -NumericalGenotypePlugin -endPlugin -export $workdir/3a_genos.numeric.txt -exportType ReferenceProbability
# 
# # Extract genotype coordinates for GAPIT map file
# echo -e "SNP\tChromosome\tPosition" > $workdir/3a_genos.map.txt
# zcat $genos | tail -n +2 | cut -f1,3,4  >>  $workdir/3a_genos.map.txt

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
    Rscript $scriptdir/3b_ResampleGapit.r --genofile $genofile --phenofile $phenofile --map $mapfile --trait $trait --num-pcs $num_pcs \
        --permdir $mypermdir --resampledir $myresampledir --num-perms $num_perms --num-resamples $num_resamples --resample-fraction $resample_fraction

#     # Save actual permutation results with p-values; remove all other intermediate files
#     pigz $mypermdir/*/GAPIT*.GWAS.Results.csv
#     rm $mypermdir/*/*.pdf $mypermdir/*/*.csv
#     
#     # Save actual resample results with p-values; remove all other intermediate files
#     pigz $myresampledir/*/GAPIT*.GWAS.Results.csv
#     rm $myresampledir/*/*.pdf $myresampledir/*/*.csv
#     
#     
#     # Collect permuted p-values
#     myperms=$permdir/$trait/*/*.GWAS.Results.csv.gz
#     myresamples=$resampledir/$trait/*/*.GWAS.Results.csv.gz   
#     Rscript 1e_GetRmips.r --perms $myperms --resamples $myresamples --perm-p-cutoff 0.05 --rmip-cutoff 0.05 -o $resultdir/1e_$trait.p05.rmip05 -t $trait
    
    break
done


# TODO: Get linkage among all the RMIP hits and see if any should be combined into single hits?
