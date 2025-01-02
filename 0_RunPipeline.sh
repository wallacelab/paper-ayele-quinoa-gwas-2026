#! /bin/bash

# Run genotyping of Quinoa SNPs

# Note: restriction enzymes were supposedly HaeIII (GG|CC) + MspI (C|CGG). HaeIII is a blunt-end cutter and all sequences start with TAA, so not sure what's going on.
#       Looking at an example BAM alignment (1001), it looks like it was a TTAA cutter, and MseI cuts T|TAA, which would match the output
#          Fwd (relative to genome)    Rev (relative to genome)
#              AGTTTAACGAT                     GGTTAAGGAT
#             GAATTTAAAAAC                     TTTTAATAA
#              TATTTAAGG                       GCTTAAAGC
#             AATTTTAAGA                      CGATTAAACT
#             AATTTTAAGG                       CATTAAAGT
#              CGATTAAAA                        CTTAATC



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



# # Basic quality control (fastqc, multiqc, adaptor trimming, etc. Looking at one of CD Genomics' BAM files, everything was still 144 bp, so no adapter trimming?)
# # NOTE: Not trimming off after a second cut site; will rely on chimera detection to filter those out
# adapters=$scriptdir/adapters.fa
# echo "" | cat $adapterdir/TruSeq3-PE-2.fa - $adapterdir/NexteraPE-PE.fa > $adapters  # Echo part is to put a newline between the files; using both to be conservative. Nextera transposon seqs definitely present.
# min_overlap=15  # Minimum overlap for joining paired ends reads with PEAR
# max_pval=0.001  # MAximum p-value for paired end joining with PEAR
# bash $scriptdir/1_QualityControl.sh $scriptdir $qcdir $rawdir $trimdir $pairdir $adapters $min_overlap $max_pval $num_cores


# Call SNPs; supposedly mpileup only uses primary alignments, so that plus min quality 20 should remove chimeras
max_depth=10000 # Maximum depth to pileup
min_base_quality=20 # Minimum base quality for pileup
subsample_for_stats=100000 # Subsample the raw SNP file to this amount to get genotyping stats
minDP=5  # Minimum genotype depth to keep a genotype call
maxDP=200 # Maximum genotype depth to keep a genotype call
min_maf=0.05   # Minimum minor allele freq to keep a site
max_missing=0.5 # Max fraction missing to keep a site
max_het=0.25 # Max p% heterozygous to keep a site
bash $scriptdir/2_CallSnps.sh $scriptdir $snpdir $trimdir $joindir $genome $max_depth $min_base_quality $subsample_for_stats \
    $minDP $maxDP $min_maf $max_missing $max_het $num_cores

