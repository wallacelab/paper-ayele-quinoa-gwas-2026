#! /bin/bash

# Run genotyping of Quinoa SNPs

genome=/home/jgwall/Projects/0_RawData/Genomes/QuinoaGenome/Cquinoa_QQ74_v2.fasta.gz
srcdir=/home/jgwall/Projects/0_RawData/Sequencing/2024_12_13_Addisu_quinoa_GBS/01.rawdata/
rawdir=0_raw_data
num_cores=7  # Number of CPU cores to use


# Set up directories for major steps
scriptdir=0_Scripts
qcdir=1_QualityControl
trimdir=$qcdir/1b_Trimmed
snpdir=2_CallSnps
if [[ ! -e $rawdir ]]; then ln -s $srcdir $rawdir; fi  # Link raw data instead of copying
if [[ ! -e $qcdir ]] ; then mkdir $qcdir; fi
if [[ ! -e $snpdir ]] ; then mkdir $snpdir; fi


# TODO: Set up basic quality control with fastq, multiqc, and other checks. Still need to make adapters file, too?

# Basic quality control (fastqc, multiqc, etc.)
adapters=$scriptdir/adapters.fa
#echo "" | cat $adapterdir/TruSeq3-PE-2.fa - $adapterdir/NexteraPE-PE.fa > $adapters  # Echo part is to put a newline between the files; using both to be conservative. Nextera transposon seqs definitely present.
bash $scriptdir/1_QualityControl.sh $scriptdir $qcdir $rawdir $trimdir $adapters $num_cores




# # Call SNPs
# bash $scriptdir/2_CallSnps.sh

