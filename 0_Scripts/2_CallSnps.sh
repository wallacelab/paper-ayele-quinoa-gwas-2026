#! /bin/bash

# Align reads to Quinoa genome and call SNPs, including quality control

# Arguments
scriptdir=$1 # Folder with accessory scripts in it
workdir=$2  # Working directory for this script
trimdir=$3  # Directory where trimmed fastq sequences need to end up
joindir=$4  # Directory where pair-joined fastq files need to end up
genome=$5   # Genome to align against (=fasta file), already indexed by BWA
max_depth=$6 # Maximum depth to pileup
min_base_quality=$7 # Minimum base quality for pileup
subsample_for_stats=$8 # Subsample the raw SNP file to this amount to get genotyping stats
minDP=$9
maxDP=${10}
min_maf=${11}   # Minimum minor allele freq to keep a site
max_missing=${12} # Max fraction missing to keep a site
max_het=${13} # Max p% heterozygous to keep a site
max_mean_depth=${14} # Max mean depth per site (to help remove pileups from paralogs)
highqual_max_missing_site=${15} # Max missing (sites) for high quality filters
highqual_max_missing_taxon=${16} # Max missing (taxon) for high quality filters
highqual_min_mean_depth=${17} # Minimum mean readdepth for high quality filters
num_cores=${18}

# Relevant folders
bamdir=$workdir/2a_alignments
snpdir=$workdir/2b_raw_snps

# Set up folder structures
if [ ! -e $workdir ]; then mkdir $workdir; fi
if [ ! -e $bamdir ]; then mkdir $bamdir; fi
if [ ! -e $snpdir ]; then mkdir $snpdir; fi

 
###################
# CONDA ENVIRONMENTS
###################

# Conda environment names. (These are handled with Miniconda and are meant to make things reproducible by making identical environments)
conda_bwa=bwa-0.7.17
conda_bcftools=bcftools-1.9
conda_samtools=samtools-1.16.1
conda_tassel5=tassel-5.2.89
conda_vcftools=vcftools-0.1.16
TASSEL5="run_pipeline.pl -Xms10g -Xmx40g" # shortcut for running TASSEL

# Load functions required to be able to activate conda environments within scripts.
. $(conda info --root)/etc/profile.d/conda.sh   # Loads the functions required to activate conda; KEEP THIS COMMAND UNCOMMENTED

###################
# Genome Alignment
##################

# # Align to genome
# for fwd in $trimdir/*_1P.fastq.gz; do
# 
#     # Derive names and files 
#     sample=`basename $fwd`
#     sample=${sample/_1P.fastq.gz/}
#     rev=${fwd/_1P.fastq.gz/_2P.fastq.gz}
# 
#     # Align
#     conda activate $conda_bwa
#     bwa mem -t $num_cores $genome $fwd $rev > $bamdir/$sample.sam     
# 
#     # Compress, sort, & index
#     conda activate $conda_samtools
#     samtools sort -@ $num_cores -o $bamdir/$sample.bam $bamdir/$sample.sam
#     samtools index $bamdir/$sample.bam
#   
#     # Cleanup
#     rm $bamdir/$sample.sam
# 
#     # break
# done



###################
# Run Genotyping
###################

# # # Note: The commands that CD Genomics originally used to call SNPs (from their VCF file) was:
# # #    samtoolsCommand=samtools mpileup -q 1 -C 50 -t SP -t DP -m 2 -F 0.002 /
# # #            -uf /TJPROJ6/CCX/DataBase/Chenopodium_quinoa/quinoa.kazusa.or.jp/Cqu_r1.0.fasta /
# # #            -b /TJNAS01/CCX/Project/Variation/X101SC19121988-Z01_an100biqi_jy/F019_part2_merge_part1/BWA/sort.bam.list  /
# # #            -l /TJNAS01/CCX/Project/Variation/X101SC19121988-Z01_an100biqi_jy/F019_part2_merge_part1/SNP/01.split/Chenopodium_quinoa_split_0
# # #    reference=file:///TJPROJ6/CCX/DataBase/Chenopodium_quinoa/quinoa.kazusa.or.jp/Cqu_r1.0.fasta

# conda activate $conda_bcftools


# bcftools mpileup is still a single-core process, apparently, so split by chromosome to be able to do parallel processing

# # Get chromosomes & contigs from genome index
# chroms=`cut -f1 $genome.fai | grep "^Cq"`
# contigs=`cut -f1 $genome.fai | grep -v "^Cq" | tr '\n' ',' | sed -r -e "s/,$//"` # Join contigs with commas and remove terminal comma
# 
# ## Make commands
# commands=$workdir/2b_commands.txt
# echo "" > $commands # Clear out from previous attempts
# for chrom in $contigs $chroms ; do
#     # determine output name
#     if [[ "$chrom" == *"Contig"* ]] ; then
#         name="contigs"
#     else
#         name=$chrom
#     fi
#     
#     # Record command for later parallel processing
#     echo "bcftools mpileup -d $max_depth -f "$genome" -r $chrom -Q $min_base_quality --output-type u --skip-indels --annotate FORMAT/AD,FORMAT/DP $bamdir/*.bam | " \
#         "bcftools call --multiallelic-caller --skip-variants indels --output $snpdir/2b_$name.raw.bcf --output-type b --variants-only " >> $commands
# done
# 
# # Run in parallel
# cat $commands | parallel --progress


# # Concatenate into a single file, chromosomes before contigs (like in genome)
# bcftools concat -o $workdir/2b_snps.raw.bcf --output-type b --threads $num_cores $snpdir/2b_Cq*.raw.bcf $snpdir/2b_contigs.raw.bcf 
# bcftools index $workdir/2b_snps.raw.bcf


###########
# Filtering
###########

conda activate $conda_bcftools

# # # # # # Count how many sites in total (can get this from index?) - DEPRECATED
# # # # # bcftools index --stats $workdir/2b_snps.raw.bcf > $workdir/2b_snps.raw.stats.txt
# # # # # Rscript -e "counts=read.delim('$workdir/2b_snps.raw.stats.txt')[,3]; cat('Total',sum(counts),'SNPs in raw genotype file')"

# # # # # # Take a random sample of X sites to run stats on - DEPRECATED
# # # # # # bcftools view $workdir/2b_snps.raw.bcf | grep -v "^#" | cut -f 1-2 | gzip > $workdir/2b_snps.raw.sites.txt.gz
# # # # # # Rscript -e "x=read.delim('$workdir/2b_snps.raw.sites.txt.gz'); set.seed(1); tokeep = sort(sample(1:nrow(x), replace=FALSE, size=$subsample_for_stats))" \
# # # # # #     -e "write.table(x[tokeep,], file='$workdir/2c_snps.raw.subset.sites.txt', row.names=FALSE, col.names=FALSE, sep='\t', quote=FALSE)"  # Random subsample and write out new sites
# # # # # # bcftools filter --regions-file $workdir/2c_snps.raw.subset.sites.txt --output-type z -o $workdir/2c_snps.raw.subset.vcf.gz $workdir/2b_snps.raw.bcf
# # # # # # bcftools index $workdir/2c_snps.raw.subset.vcf.gz

# # Convert to VCF so TASSEL & vcftools can read it
# bcftools view  --output-type z -o $workdir/2c_snps.raw.vcf.gz $workdir/2b_snps.raw.bcf

# Function to generate all the stats and data for plotting diagnostics; arguments are input gzipped VCF file and output plot PNG
# NOTE: Resuting plot should be used to adjust the filter parameters in the master 0_RunPipeline.sh script
plot_diagnostics () {
    # Define variables
    input=$1
    statdir=$1.data
    outfile=$2
    if [ ! -e $statdir ]; then mkdir $statdir; fi
    
    # TASSEL genotype and taxa summaries
    conda activate $conda_tassel5
    $TASSEL5 -SortGenotypeFilePlugin -inputFile $input -outputFile $statdir/sorted.vcf.gz
    $TASSEL5 -vcf $statdir/sorted.vcf.gz -genotypeSummary site -export $statdir/sitesummary.txt
    $TASSEL5 -vcf $statdir/sorted.vcf.gz -genotypeSummary taxa -export $statdir/taxasummary.txt
    
    # VCF tools depth summaries
    conda activate $conda_vcftools
    vcftools --gzvcf $statdir/sorted.vcf.gz  --site-mean-depth --stdout > $statdir/site_depth.txt
    vcftools --gzvcf $statdir/sorted.vcf.gz  --geno-depth --stdout | tail -n+2 | cut -f 3- | pigz > $statdir/geno_depth.txt.gz  # skip header and 1st 2 columns so is just depth info  
    
    # Plot diagnostics
    Rscript $scriptdir/2c_PlotGenoSummary.r --depthfile $statdir/site_depth.txt --genodepth $statdir/geno_depth.txt.gz \
        --sitefile $statdir/sitesummary.txt --taxafile $statdir/taxasummary.txt \
        --max-datapoints 10000 -o $outfile
}


# # Filter step 1 - Genotype depth (but also getting other stats to compare)
# plot_diagnostics $workdir/2c_snps.raw.vcf.gz $workdir/2c_summary_plots.raw.png
# 
# # Filter genotypes based on genotype depth first
# conda activate $conda_vcftools
# vcftools --gzvcf $workdir/2c_snps.raw.vcf.gz  --minDP $minDP --maxDP $maxDP --recode --stdout | pigz > $workdir/2d_snps.genodepth_${minDP}_${maxDP}.vcf.gz
# plot_diagnostics $workdir/2d_snps.genodepth_${minDP}_${maxDP}.vcf.gz $workdir/2d_summary_plots.genodepth_${minDP}_${maxDP}.png 
 

# # Filter step 2 - MAF, missingness, hets, etc.
# # Filter based on diagnostic plots
# statdir=$workdir/2d_snps.genodepth_${minDP}_${maxDP}.vcf.gz.data
# Rscript $scriptdir/2d_GetFilterLists.r --sitefile $statdir/sitesummary.txt --taxafile $statdir/taxasummary.txt --depthfile $statdir/site_depth.txt \
#     --site-min-maf $min_maf --site-max-missing $max_missing --site-max-het $max_het --site-max-mean-depth $max_mean_depth \
#     --outtaxa $workdir/2d_taxa_to_keep.txt --outsites $workdir/2d_sites_to_keep.txt
# 
# # TASSEL apparently capitalizes the chromosome names, so need to fix that
# cat $workdir/2d_sites_to_keep.txt | sed -r -e "s|^CQ|Cq|" -e "s|^CONTIG|Contig|" > $workdir/2d_sites_to_keep.fix_case.txt
#  
# # conda activate $conda_bcftools
# bcftools view --targets-file $workdir/2d_sites_to_keep.fix_case.txt --samples-file $workdir/2d_taxa_to_keep.txt --output-type z --threads $num_cores \
#     --output-file $workdir/2e_snps.filtered.vcf.gz $workdir/2d_snps.genodepth_${minDP}_${maxDP}.vcf.gz
# plot_diagnostics $workdir/2e_snps.filtered.vcf.gz $workdir/2e_summary_plots.filtered.png


# # Filter step 3 - Higher quality (= more stringent) filtering. Not sure if will need
# statdir=$workdir/2e_snps.filtered.vcf.gz.data
# Rscript $scriptdir/2d_GetFilterLists.r --sitefile $statdir/sitesummary.txt --taxafile $statdir/taxasummary.txt --depthfile $statdir/site_depth.txt \
#     --site-max-missing $highqual_max_missing_site --taxa-max-missing $highqual_max_missing_taxon --site-min-mean-depth $highqual_min_mean_depth \
#     --outtaxa $workdir/2f_taxa_to_keep.highqual.txt --outsites $workdir/2f_sites_to_keep.highqual.txt
# 
# # TASSEL apparently capitalizes the chromosome names, so need to fix that
# cat $workdir/2f_sites_to_keep.highqual.txt | sed -r -e "s|^CQ|Cq|" -e "s|^CONTIG|Contig|" > $workdir/2f_sites_to_keep.highqual.fix_case.txt
#  
# # conda activate $conda_bcftools
# bcftools view --targets-file $workdir/2f_sites_to_keep.highqual.fix_case.txt --samples-file $workdir/2f_taxa_to_keep.highqual.txt --output-type z --threads $num_cores \
#     --output-file $workdir/2f_snps.filtered.highqual.vcf.gz $workdir/2e_snps.filtered.vcf.gz
plot_diagnostics $workdir/2f_snps.filtered.highqual.vcf.gz $workdir/2f_summary_plots.filtered.highqual.png




# NOTE The newest Quinoa genome (QQ74_v2, from Rey et al 2023) did a diversity study where "only biallelic variants with minimum read depth 5, maximum missing data 20% and MAF 0.01% were retained"




# # # # 
# # # # 
# # # # # ##############
# # # # # # Cladogram
# # # # # ##############
# # # # # 
# # # # # # Tree
# # # # # $TASSEL5 -vcf 1g_hemp_filtered.sorted.vcf.gz -tree Neighbor -treeSaveDistance false -export 1g_hemp_filtered.sorted.tre -exportType Text
# # # # # mv 1g_hemp_filtered.sorted.tre.txt 1g_hemp_filtered.sorted.tre
# # # # 
# # # # 
