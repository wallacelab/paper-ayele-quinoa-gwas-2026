#! /bin/bash

 
# Note: The commands that were originally used to call SNPs (from their VCF file) was:
#    samtoolsCommand=samtools mpileup -q 1 -C 50 -t SP -t DP -m 2 -F 0.002 /
#            -uf /TJPROJ6/CCX/DataBase/Chenopodium_quinoa/quinoa.kazusa.or.jp/Cqu_r1.0.fasta /
#            -b /TJNAS01/CCX/Project/Variation/X101SC19121988-Z01_an100biqi_jy/F019_part2_merge_part1/BWA/sort.bam.list  /
#            -l /TJNAS01/CCX/Project/Variation/X101SC19121988-Z01_an100biqi_jy/F019_part2_merge_part1/SNP/01.split/Chenopodium_quinoa_split_0
#    reference=file:///TJPROJ6/CCX/DataBase/Chenopodium_quinoa/quinoa.kazusa.or.jp/Cqu_r1.0.fasta



###################
# CREATE CONDA ENVIRONMENT
###################

# Conda environment names. (These are handled with Miniconda and are meant to make things reproducible by making identical environments)
conda_bcftools=bcftools-1.9

# Load functions required to be able to activate conda environments within scripts.
. $(conda info --root)/etc/profile.d/conda.sh   # Loads the functions required to activate conda; KEEP THIS COMMAND UNCOMMENTED
conda activate $conda_bcftools


####################
# Quality Control
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



###################
# Run Genotyping
###################

# # Index genome for alignment (bwa version 0.7.12-r1039)
# bwa index $genome 

# # Align each sample to the genome
# for sample in $src/*.fastq.gz; do
#     name=`basename $sample`
#     name=${name/_R1*/}
#     echo "Processing alignments for $name"
#     
# #     Align to genome
#     bwa aln -t $numcores -f $workdir/1_$name.sai $genome $sample
# 
#     # Convert to SAM
#     bwa samse $genome $workdir/1_$name.sai $sample > $workdir/1_$name.sam
#     
# #     Convert to BAM & sort
#     samtools view -h -S -b $workdir/1_$name.sam > $workdir/1_$name.bam
#     samtools sort -@ $numcores $workdir/1_$name.bam $workdir/1_$name.sorted
#     samtools index $workdir/1_$name.sorted.bam
# 
#     rm $workdir/1_$name.sam $workdir/1_$name.sai $workdir/1_$name.bam 
# #     break
# done

# # Index genome for samtools
# samtools faidx $genome  

# # Get list of locations to call at; have to tweak the format
# zcat $public_snps | vcf2bed | cut -f1-3 > $workdir/1a_target_locations.bed
# cat $workdir/1a_target_locations.bed | sed -r -e "s/^/Chr0/" | sed -e "s/Chr010/Chr10/" > $workdir/1a_target_locations.corrected.bed

# # Call SNPs on all samples (skip indels with --skip-variants b/c there are no indels in the Hu et al data and it causes issues trying to merge & impute them)
# max_depth=1000
# min_base_quality=20 
# bcftools mpileup -d $max_depth -f $genome --regions-file $workdir/1a_target_locations.corrected.bed -Q $min_base_quality --output-type z --output $workdir/1b_raw_probs.vcf.gz  $workdir/*.sorted.bam  
# bcftools call --multiallelic-caller --threads $numcores --skip-variants indels --output $workdir/1c_snp_calls.vcf.gz --output-type z $workdir/1b_raw_probs.vcf.gz 

# # Replace sample, site, and chromosome names throughout the VCF file.
# python3 1d_RenameVcfSamplesAndSites.py -i $workdir/1c_snp_calls.vcf.gz -o $workdir/1d_snp_calls.renamed.vcf.gz --sample-keyfile $sample_name_key # --debug


# # Make keys to correct sample and site names for when merged with existing key

# DEPRECATED zcat $workdir/1c_snp_calls.vcf.gz | grep "^Chr" | sed -r "s|^Chr(..)\t([0-9]+)\t.+|Chr\1\t\2\tS\1_\2|" | sed "s|\tS0|\tS|"  > $workdir/1d.site-rename-key.txt  # Capture chrom + position to make a site name compatible with existing data




# Details from Hu et al for how they got the alignment
# Do I need to do GBS, or just realign and call at these specific locations? I think the latter makes more sense.
# Burrows–Wheeler Alignment version 0.7.10 (Li and Durbin, 2009). The raw GBS
# reads were converted to tags with the -GBSSeqToTag-DBPlugin tool with -kmerLength 64, -minKmerL 20,
# -mnQS 20, and -mxKmerNum 100000000. The tags were converted to fastq format with -TagExportToFastqPlugin,
# then aligned to the reference genome using aln option and converted to SAM format using samse option in
# Burrows–Wheeler Alignment. 
