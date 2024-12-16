#! /bin/bash

# Do basic QC on sequences
# Note: Trimmomatic parameters were done to match KBase defaults, except for adding parameters to remove adapters.
#       (I Spot-checked results before adding apadter-trim and confirmed were the same)

# Arguments
scriptdir=$1 # Folder with accessory scripts in it
workdir=$2  # Working directory for this script
datadir=$3  # Raw data with original sequencing files
trimdir=$4  # Directory where trimmed fastq sequences need to end up
adapters=$5 # Sequencing adapters in a FASTA file
num_cores=$6

# Relevant folders
fastqcdir=$workdir/1a_FastQC

# Set up folder structures
if [ ! -e $workdir ]; then mkdir $workdir; fi
if [ ! -e $fastqcdir ]; then mkdir $fastqcdir; fi
if [ ! -e $trimdir ]; then mkdir $trimdir; fi
if [ ! -e $trimdir.fastqc ]; then mkdir $trimdir.fastqc; fi


###################
# CONDA ENVIRONMENTS
###################

# Conda environment names. (These are handled with Miniconda and are meant to make things reproducible by making identical environments)
conda_fastqc=fastqc-0.11.9      # Conda environment with FASTQC and MultiQC
conda_trimmomatic=trimmomatic-0.39

# Load functions required to be able to activate conda environments within scripts.
. $(conda info --root)/etc/profile.d/conda.sh   # Loads the functions required to activate conda; KEEP THIS COMMAND UNCOMMENTED

##############
# FASTQC of raw data
##############

# FastQC analysis of raw samples
conda activate $conda_fastqc
fastqc -f fastq -o $fastqcdir -t $num_cores $datadir/*.fq.gz
multiqc $fastqcdir -o $fastqcdir.multiqc


# ##############
# # Adapter trimming and paired end joining
# ##############
# 
# # Cut adapters
# conda activate $conda_trimmomatic
# for fwd in $datadir/*_R1_001.fastq.gz; do
#     rev=${fwd/_R1_001/_R2_001}
#     sample=${fwd/_R1_001.fastq.gz/}
#     sample=${sample/$datadir\//}
#     trimmomatic PE -threads $num_cores -basein $fwd -baseout $trimdir/$sample.fastq.gz \
#         ILLUMINACLIP:$adapters:2:30:10:5:True SLIDINGWINDOW:4:15 LEADING:3 TRAILING:3 MINLEN:36 TOPHRED33
# #     break
# done

# ##############
# # FASTQC of trimmed data
# ##############
# 
# # FastQC analysis of trimmed data
# conda activate $conda_fastqc
# fastqc -f fastq -o $trimdir.fastqc -t $num_cores $trimdir/*P.fastq.gz  # *P is specifically the paired ones as opposed to unpaired ones
# multiqc $trimdir.fastqc -o $trimdir.multiqc



 
# Note: The commands that were originally used to call SNPs (from their VCF file) was:
#    samtoolsCommand=samtools mpileup -q 1 -C 50 -t SP -t DP -m 2 -F 0.002 /
#            -uf /TJPROJ6/CCX/DataBase/Chenopodium_quinoa/quinoa.kazusa.or.jp/Cqu_r1.0.fasta /
#            -b /TJNAS01/CCX/Project/Variation/X101SC19121988-Z01_an100biqi_jy/F019_part2_merge_part1/BWA/sort.bam.list  /
#            -l /TJNAS01/CCX/Project/Variation/X101SC19121988-Z01_an100biqi_jy/F019_part2_merge_part1/SNP/01.split/Chenopodium_quinoa_split_0
#    reference=file:///TJPROJ6/CCX/DataBase/Chenopodium_quinoa/quinoa.kazusa.or.jp/Cqu_r1.0.fasta

# The newest Quinoa genome (QQ74_v2, from Rey et al 2023) did a diversity study where "only biallelic variants with minimum read depth 5, maximum missing data 20% and MAF 0.01% were retained"

# ###################
# # CREATE CONDA ENVIRONMENT
# ###################
# 
# # Conda environment names. (These are handled with Miniconda and are meant to make things reproducible by making identical environments)
# conda_bcftools=bcftools-1.9
# 
# # Load functions required to be able to activate conda environments within scripts.
# . $(conda info --root)/etc/profile.d/conda.sh   # Loads the functions required to activate conda; KEEP THIS COMMAND UNCOMMENTED
# conda activate $conda_bcftools
# 
# 
# ####################
# # Quality Control
####################




# # # # ####################
# # # # # Check 2 - Make sure sequencing is actually sorghum. (Initial alignments only had ~20% alignment with BWA
# # # # ####################
# # # # 
# # # # 
# # # # # Software
# # # # KRAKEN2="$HOME/Software/Metagenomics/kraken2/bin/kraken2"    # Kraken 2.0.7-beta
# # # # KRONA="$HOME/Software/Metagenomics/KronaTools-2.7/bin/ktImportText" # KronaTools 2.7
# # # # 
# # # # # Files
# # # # datadir=$HOME/Projects/0_RawData/Sequencing/2019_11_20_Sorghum_GBS_for_FVSU/
# # # # kraken_db="/home/jgwall/Software/Metagenomics/kraken2/kraken_plants"
# # # # 
# # # # # Working directory
# # # # workdir=0b_KrakenCheck
# # # # if [ ! -e $workdir ]; then mkdir $workdir; fi
# # # # 
# # # # maxprocs=7 # Max number of threads to use
# # # # for infile in $datadir/*.fastq.gz; do
# # # #     name=`basename $infile`
# # # #     name=${name/_R1*/}
# # # # 
# # # #     # Align to Kraken database 
# # # #     $KRAKEN2 --db $kraken_db --threads $maxprocs --report $workdir/0b_$name.kraken_results.summary.txt --use-names $infile
# # # #     
# # # #     #break
# # # # done
# # # # 
# # # # # Everything is majority Sorghum, ranging fro 80% to 95% (with most of the rest unclassified), so that looks good.



# #! /bin/bash
# 
# # Get SNPs for each of the samples
# 
# # Add BedOps package to file path
# BEDOPS=/home/jgwall/Software/FileTools/BedOps
# export PATH="$PATH:$BEDOPS" 
# 
# genome=$HOME/Projects/0_RawData/SorghumGenome/3.1_phytozome/Sbicolor_313_v3.0.fa
# src=$HOME/Projects/0_RawData/Sequencing/2019_11_20_Sorghum_GBS_for_FVSU/
# public_snps=$HOME/Projects/0_RawData/SorghumGermplasm/Hu_et_al_10k_GBS_genotypes/SNPs_imp.recode.vcf.bgzip.gz
# sample_name_key=0_sample_rename_key.txt
# numcores=7 # Number of cores/threads to use
# 
# workdir=1_CallGenotypes
# if [ ! -e $workdir ]; then mkdir $workdir; fi


