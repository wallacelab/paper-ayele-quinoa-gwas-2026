#! /bin/bash

# Run genotyping of Quinoa SNPs

# Note: restriction enzymes were HaeIII (GG|CC) + MspI (C|CGG). HaeIII is a blunt-end cutter and all sequences start with TAA, so not sure what's going on.


# Set up global variables
genome=/home/jgwall/Projects/0_RawData/Genomes/QuinoaGenome/Cquinoa_QQ74_v2.fasta.gz  # Genome
srcdir=/home/jgwall/Projects/0_RawData/Sequencing/2024_12_13_Addisu_quinoa_GBS/01.rawdata/
adapterdir="/home/jgwall/Software/Python/miniconda3/envs/trimmomatic-0.39/share/trimmomatic-0.39-2/adapters/" # Location of trimmomatic adapters
rawdir=0_raw_data
num_cores=7  # Number of CPU cores to use


# Set up directories for major steps
scriptdir=0_Scripts
qcdir=1_QualityControl
trimdir=$qcdir/1b_Trimmed
joindir=$qcdir/1c_PairJoining
snpdir=2_CallSnps
if [[ ! -e $rawdir ]]; then ln -s $srcdir $rawdir; fi  # Link raw data instead of copying
if [[ ! -e $qcdir ]] ; then mkdir $qcdir; fi
if [[ ! -e $snpdir ]] ; then mkdir $snpdir; fi


# TODO - Remove read merging. samtools mpileup handles overlaps automatically


# # Basic quality control (fastqc, multiqc, adaptor trimming, etc. Looking at one of CD Genomics' BAM files, everything was still 144 bp, so no adapter trimming?)
# # NOTE: Not trimming off after a second cut site; will rely on chimera detection to filter those out
# adapters=$scriptdir/adapters.fa
# echo "" | cat $adapterdir/TruSeq3-PE-2.fa - $adapterdir/NexteraPE-PE.fa > $adapters  # Echo part is to put a newline between the files; using both to be conservative. Nextera transposon seqs definitely present.
# min_overlap=15  # Minimum overlap for joining paired ends reads with PEAR
# max_pval=0.001  # MAximum p-value for paired end joining with PEAR
# bash $scriptdir/1_QualityControl.sh $scriptdir $qcdir $rawdir $trimdir $pairdir $adapters $min_overlap $max_pval $num_cores

# TODO: When align, find way of removing chimeras. BWA can report 2+ alignments for chimeras, so maybe just choose singly mapping ones?


# Call SNPs
bash $scriptdir/2_CallSnps.sh $scriptdir $snpdir $trimdir $joindir $genome $num_cores

