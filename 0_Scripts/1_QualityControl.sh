#! /bin/bash

# Do basic QC on sequences
# Note: Trimmomatic parameters were done to match KBase defaults, except for adding parameters to remove adapters.
#       (I Spot-checked results before adding apadter-trim and confirmed were the same)

# Arguments
scriptdir=$1 # Folder with accessory scripts in it
workdir=$2  # Working directory for this script
datadir=$3  # Raw data with original sequencing files
trimdir=$4  # Directory where trimmed fastq sequences need to end up
pairdir=$5  # Directory where pair-joined fastq files need to end up
adapters=$6 # Sequencing adapters in a FASTA file
min_overlap=$7 # Minimum overlap for joining paired ends reads with PEAR
max_pval=$8  # MAximum p-value for paired end joining with PEAR
num_cores=$9

# Relevant folders
fastqcdir=$workdir/1a_FastQC

# Set up folder structures
if [ ! -e $workdir ]; then mkdir $workdir; fi
if [ ! -e $fastqcdir ]; then mkdir $fastqcdir; fi
if [ ! -e $trimdir ]; then mkdir $trimdir; fi
if [ ! -e $trimdir.fastqc ]; then mkdir $trimdir.fastqc; fi
if [ ! -e $pairdir ]; then mkdir $pairdir; fi


###################
# CONDA ENVIRONMENTS
###################

# Conda environment names. (These are handled with Miniconda and are meant to make things reproducible by making identical environments)
conda_fastqc=fastqc-0.11.9         # Conda environment with FASTQC and MultiQC
conda_trimmomatic=trimmomatic-0.39 # Trimmomatic
conda_pear=pear-0.9.6              # PEAR paired-end read merger

# Load functions required to be able to activate conda environments within scripts.
. $(conda info --root)/etc/profile.d/conda.sh   # Loads the functions required to activate conda; KEEP THIS COMMAND UNCOMMENTED

# ##############
# # FASTQC of raw data
# ##############
# 
# # FastQC analysis of raw samples
# conda activate $conda_fastqc
# fastqc -f fastq -o $fastqcdir -t $num_cores $datadir/*.fq.gz
# multiqc $fastqcdir -o $fastqcdir.multiqc


# ##############
# # Adapter trimming
# ##############
# 
# # Cut adapters
# conda activate $conda_trimmomatic
# for fwd in $datadir/*_1.fq.gz; do
#     rev=${fwd/_1.fq.gz/_2.fq.gz}
#     sample=${fwd/_1.fq.gz/}
#     sample=${sample/$datadir\//}
#     trimmomatic PE -threads $num_cores -basein $fwd -baseout $trimdir/$sample.fastq.gz \
#         ILLUMINACLIP:$adapters:2:30:10:5:True SLIDINGWINDOW:4:20 LEADING:3 TRAILING:3 MINLEN:36
#     #break
# done

# ##############
# # FASTQC of trimmed data
# ##############
# 
# # FastQC analysis of trimmed data
# conda activate $conda_fastqc
# fastqc -f fastq -o $trimdir.fastqc -t $num_cores $trimdir/*P.fastq.gz  # *P is specifically the paired ones as opposed to unpaired ones
# multiqc $trimdir.fastqc -o $trimdir.multiqc

# # Plot read retention
# Rscript $scriptdir/1b_PlotReadRetention.r -a $fastqcdir.multiqc/multiqc_data/multiqc_general_stats.txt -b $trimdir.multiqc/multiqc_data/multiqc_general_stats.txt \
#     -o $workdir/1b_read_retention_after_trimming


###################
# Join paired ends - DEPRECATED. Samtools mpileup will remove overlaps by default, so no need to merge
###################

# # Join reads
# conda activate $conda_pear
# for fwd in $trimdir/*_1P.fastq.gz; do
#     rev=${fwd/_1P/_2P}
#     sample=${fwd/_1P.fastq.gz/}
#     sample=${sample/$trimdir\//}
#     pear  -f $fwd -r $rev -o $pairdir/$sample --min-overlap $min_overlap --p-value $max_pval --threads $num_cores
#     pigz $pairdir/$sample*.fastq # Zip up output files
# done

# # Count reads in each class
# countfile=$workdir/1c_pairjoin_file_stats.txt
# echo -e "sample\tset\treads" > $countfile
# for fastq in $pairdir/*.fastq.gz; do
#     # File info
#     info=`basename $fastq`
#     info=${info/.fastq.gz/}
#     sample=${info/\.*/}
#     set=${info/*\./}
#     
#     # Read counts
#     count=`zcat $fastq | wc -l`
#     count=$(( $count / 4 ))
#     
#     echo -e "$sample\t$set\t$count" >> $countfile
#     #break
# done
# 
# Rscript $scriptdir/1c_PlotReadJoinStats.r -i $countfile -o $workdir/1c_read_join_stats

 
