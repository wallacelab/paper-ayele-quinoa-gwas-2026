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
perm_p_cutoff="$9"  # P-value cutoff for permutations to count as null "hits"
rmip_cutoff="${10}" # RMIP cutoff for a hit to count as real
genome=${11}
models="${12}" # Which models to run in GAPIT (BLINK, FarmCPU, etc)
num_perms_glm=${13}  # Number of permutations for GLM empirical p-values
num_cores=${14} # Number of CPU cores to use

# Manhattan plot parameters (not important enough to pass from master script)
winsize=5e6
step=5e5

# Make output directory if don't have yet
tasseldir=$workdir/3g_TASSEL_GLM
if [[ ! -e $workdir ]] ; then mkdir $workdir; fi
if [[ ! -e $tasseldir ]] ; then mkdir $tasseldir; fi

########
# Conda envs
########

# Conda environments
conda_tassel5=tassel-5.2.89
conda_samtools=samtools-1.16.1

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

# # Get list of traits
# traits=`head -n 1 $workdir/3a_phenos.csv | tr ',' '\n' | tail -n +2 | tr '"' ' '`
# 
# # Key files (putting here to make clearer)
# phenofile=$workdir/3a_phenos.csv
# genofile=$workdir/3a_genos.numeric.GAPIT.rds
# mapfile=$workdir/3a_genos.map.txt
# kinfile=$workdir/3a_kinship.txt
#  
# # Loop over models
# for model in $models; do
# 
#     # Make subdirectories
#     permdir=$workdir/3b_Permutations.$model
#     resampledir=$workdir/3d_Resampling.$model
#     resultdir=$workdir/3e_Results.$model
#     if [ ! -e $permdir ] ; then mkdir $permdir; fi
#     if [ ! -e $resampledir ] ; then mkdir $resampledir; fi
#     if [ ! -e $resultdir ] ; then mkdir $resultdir; fi
# 
#     # Loop over traits
#     for trait in $traits; do
#     
#         # Output directories
#         mypermdir=$permdir/$trait
#         myresampledir=$resampledir/$trait
#         if [ ! -e $mypermdir ]; then mkdir $mypermdir; fi
#         if [ ! -e $myresampledir ]; then mkdir $myresampledir; fi
#         
#         # Run permutations & resamplings
#         Rscript $scriptdir/3b_ResampleGapit.r --genofile $genofile --phenofile $phenofile --map $mapfile --kinshipfile $kinfile --trait $trait --num-pcs $num_pcs \
#             --permdir $mypermdir --resampledir $myresampledir --num-perms $num_perms --num-resamples $num_resamples --resample-fraction $resample_fraction \
#             --model $model --num-cores $num_cores
# 
#         # Save actual permutation results with p-values; remove all other intermediate files
#         pigz $mypermdir/*/GAPIT*GWAS_Results.$model.$trait.csv
#         pigz $mypermdir/*/*Filter_GWAS_results.csv  # Probably don't need this, but keeping just in case
#         rm $mypermdir/*/*.pdf $mypermdir/*/*.csv
#         
#         # Save actual resample results with p-values; remove all other intermediate files
#         pigz $myresampledir/*/GAPIT*GWAS_Results.$model.$trait.csv
#         pigz $myresampledir/*/*Filter_GWAS_results.csv  # Probably don't need this, but keeping just in case
#         rm $myresampledir/*/*.pdf $myresampledir/*/*.csv
#          
#         
#         # Collect permuted p-values
#         for perm_p in $perm_p_cutoff; do
#             for rmip_cut in $rmip_cutoff; do
#         
#                 myperms=$mypermdir/*/GAPIT*GWAS_Results.$model.$trait.csv.gz
#                 myresamples=$myresampledir/*/GAPIT*GWAS_Results.$model.$trait.csv.gz   
#                 outprefix=$resultdir/3e_$trait.p$perm_p.rmip$rmip_cut
#                 Rscript $scriptdir/3e_GetRmips.r --perms $myperms --resamples $myresamples --perm-p-cutoff $perm_p  --rmip-cutoff $rmip_cut \
#                     -o $outprefix -t $trait --num-cores $num_cores
#             done
#         done
# 
#         #break
#     done
#     #break
# done
# 
# 
# # Make file of chromosome offsets for plotting Manhattan plots
# conda activate $conda_samtools
# zcat $genome > $workdir/tmp_genome.fa
# samtools faidx $workdir/tmp_genome.fa --fai-idx $workdir/3f_genome_index.txt
# rm $workdir/tmp_genome.fa
# 
# # Plot combined Manhattan plots
# for model in $models; do
#     for perm_p in $perm_p_cutoff; do
#         for rmip_cut in $rmip_cutoff; do
#             Rscript $scriptdir/3f_MakeManhattanPlots.r -i $workdir/3e_Results.$model/*p$perm_p.rmip$rmip_cut.csv --winsize $winsize --step $step --drop-contigs \
#                 --index $workdir/3f_genome_index.txt --rmip-cutoff $rmip_cutoff -o $workdir/3f_manahttans.$model.p$perm_p.rmip$rmip_cut
#         done
#     done
# done


###############
# GLM permutations for Manhattans with empirical p-values
###############

# conda activate $conda_tassel5
# 
# # Quick scree plot to confirm number of PCs is correct; eigenvalues are in file 2
# run_pipeline.pl -h $genos -distanceMatrix -MultiDimensionalScalingPlugin -axes 20 -endPLugin -export $tasseldir/3g_genos_mds.txt
# Rscript -e "library(ggplot2); eigens=read.delim('$tasseldir/3g_genos_mds2.txt'); myplot = ggplot(eigens) + aes(x=PC, y=eigenvalue) + geom_point()" \
#     -e "ggsave(myplot, file='$tasseldir/3g_genos_screeplot.png')"
# 
# # Take first num_pcs PCs to work with
# num_cols=$(($num_pcs + 1))
# cat $tasseldir/3g_genos_mds1.txt | cut -f 1-$num_cols | sed -r -e "s|^\S+/||g" > $tasseldir/3g_pcs.txt
# 
# # Fix genotype names so don't have path information
# zcat $genos | head -n 1 | sed -r -e "s|\t\S+/|\t|g" > $tasseldir/3g_genos.fixnames.hmp.txt
# zcat $genos | tail -n +21 >> $tasseldir/3g_genos.fixnames.hmp.txt
# pigz $tasseldir/3g_genos.fixnames.hmp.txt
# 
# # # Reformat GAPIT phenos for TASSEL
# cat $phenos | sed -r -e "s|<Taxa>|<Trait>|" > $tasseldir/3g_phenos.tassel.txt
# 
# # Split phenos for TASSEL
# ncols=`head -n 3 $tasseldir/3g_phenos.tassel.txt | tail -n 1 | tr '\t' '\n' | wc -l`
# commands=$workdir/3g_tassel_commands.txt
# echo "echo 'Starting GLM permutation analysis'" > $commands
# for traitcol in `seq 2 $ncols`; do
#     cut -f 1,$traitcol $tasseldir/3g_phenos.tassel.txt > $tasseldir/3g_phenos.tassel.$traitcol.txt
#     
#     echo "run_pipeline.pl -fork1 -h $tasseldir/3g_genos.fixnames.hmp.txt.gz -fork2 -t $tasseldir/3g_phenos.tassel.$traitcol.txt -fork3 -t $tasseldir/3g_pcs.txt \
#         -combine4 -input1 -input2 -input3 -intersect \
#         -FixedEffectLMPlugin -saveToFile false \
#         -permute true -nperm $num_perms_glm -maxP 0.1 \
#         -endPlugin \
#         -export $tasseldir/3h_glm_results_$traitcol" >> $commands
# done
# cat $commands | parallel --progress
#  
# # Rename the TASSEL results so more intuitive, then zip
# for traitcol in `seq 2 $ncols`; do
#     mv $tasseldir/3h_glm_results_${traitcol}1.txt $tasseldir/3h_glm_results_$traitcol.sites.txt
#     mv $tasseldir/3h_glm_results_${traitcol}2.txt $tasseldir/3h_glm_results_$traitcol.alleles.txt
#     
#     pigz $tasseldir/3h_glm_results_$traitcol.*.txt
# done


# Combine GLM results with RMIP results for final manhattan plots
for model in $models; do
    for perm_p in $perm_p_cutoff; do
        for rmip_cut in $rmip_cutoff; do
             Rscript $scriptdir/3i_MakeCompositeManhattanPlots.r  --rmip-infiles $workdir/3e_Results.$model/*p$perm_p.rmip$rmip_cut.csv \
                --winsize $winsize --step $step --rmip-cutoff $rmip_cutoff \
                --glm-infiles $tasseldir/3h_glm_results_*.sites.txt.gz --pval-cutoff $perm_p \
                -o $workdir/3i_combined_manahttans.$model.p$perm_p.rmip$rmip_cut 
        done
    done
done



    
