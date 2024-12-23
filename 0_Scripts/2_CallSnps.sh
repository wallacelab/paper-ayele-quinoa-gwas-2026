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
num_cores=$8

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



# bcftools mpileup is still a single-core process, apparently, so split by chromosome to be able to do parallel processing

# Get chromosomes & contigs from genome index
chroms=`cut -f1 $genome.fai | grep "^Cq"`
contigs=`cut -f1 $genome.fai | grep -v "^Cq" | tr '\n' ',' | sed -r -e "s/,$//"` # Join contigs with commas and remove terminal comma

## Make commands
commands=$workdir/2b_commands.txt
echo "" > $commands # Clear out from previous attempts
for chrom in $contigs $chroms ; do
    # determine output name
    if [[ "$chrom" == *"Contig"* ]] ; then
        name="contigs"
    else
        name=$chrom
    fi
    
    # Record command for later parallel processing
    echo "bcftools mpileup -d $max_depth -f "$genome" -r $chrom -Q $min_base_quality --output-type u --skip-indels --annotate FORMAT/AD,FORMAT/DP $bamdir/*.bam | " \
        "bcftools call --multiallelic-caller --skip-variants indels --output $snpdir/2b_$name.raw.bcf --output-type b --variants-only " >> $commands
done

# Run in parallel
conda activate $conda_bcftools
cat $commands | parallel --progress



# NOTE The newest Quinoa genome (QQ74_v2, from Rey et al 2023) did a diversity study where "only biallelic variants with minimum read depth 5, maximum missing data 20% and MAF 0.01% were retained"






# # # # ##############
# # # # # Basic Filtering - Depth
# # # # ##############
# # # # 
# # # # # bcftools view 1b_hemp.snp_calls.bcf --output-type z > 1b_hemp.snp_calls.vcf.gz
# # # # # 
# # # # # # Get reports
# # # # # $TASSEL5 -SortGenotypeFilePlugin -inputFile 1b_hemp.snp_calls.vcf.gz -outputFile 1b_hemp.snp_calls.sorted.vcf.gz
# # # # # $TASSEL5 -vcf 1b_hemp.snp_calls.sorted.vcf.gz -genotypeSummary site -export $filtdir/1c_sitesummary.txt
# # # # # $TASSEL5 -vcf 1b_hemp.snp_calls.sorted.vcf.gz -genotypeSummary taxa -export $filtdir/1c_taxasummary.txt
# # # # # vcftools --gzvcf 1b_hemp.snp_calls.sorted.vcf.gz  --site-depth --stdout > $filtdir/1c_site_depth.txt
# # # # 
# # # # # # Complicated bash one-liner to count size of predicted restriction fragments in the hemp genome
# # # # # cat "$genome" | sed -r -e "s|^>.+|=====|" | tr -d '\n' | tr [:lower:] [:upper:] | \
# # # # #     sed -r -e "s|(CA[CT]..)(..[AG]TG)|\1\n\2|g" -e "s|=====||g" | \
# # # # #     awk '{print length}' > $filtdir/1c_restriction_frag_lengths.txt
# # # # 
# # # # # # Plot diagnostic plots from the above summary data
# # # # # Rscript $scriptdir/1c_PlotGenoSummary.r --sitefile $filtdir/1c_sitesummary.txt --taxafile $filtdir/1c_taxasummary.txt --depthfile $filtdir/1c_site_depth.txt \
# # # # #     --restriction-frags $filtdir/1c_restriction_frag_lengths.txt --max-datapoints 10000 -o $filtdir/1c_summary_plots.png
# # # # 
# # # # 
# # # # # Hmm...based on number of restriction fragments and read depth, I'd expect a true mean depth of single-copy DNA around 5, but actual cluster seems around 100. 
# # # # #   (Ignoring the cluster of aligment errors that only have 1 read) Let's filter for the 100x cluster and see how it works
# # # # 
# # # # 
# # # # # # Initial filter for mean depth and % presence
# # # # # site_min_mean_depth=15
# # # # # site_max_mean_depth=125
# # # # # site_max_missing=0.2  # Equivalent to 80% present
# # # # # Rscript $scriptdir/1d_GetFilterLists.r --sitefile $filtdir/1c_sitesummary.txt --taxafile $filtdir/1c_taxasummary.txt --depthfile $filtdir/1c_site_depth.txt \
# # # # #     --site-min-mean-depth $site_min_mean_depth --site-max-mean-depth $site_max_mean_depth --site-max-missing $site_max_missing \
# # # # #     --outtaxa $filtdir/1d_taxa_to_keep.txt --outsites $filtdir/1d_sites_to_keep.txt
# # # # 
# # # # # # Perform actual filtering
# # # # # bcftools view --targets-file $filtdir/1d_sites_to_keep.txt --samples-file $filtdir/1d_taxa_to_keep.txt --output-type z --output-file $filtdir/1e_hemp_filter1.vcf.gz 1b_hemp.snp_calls.bcf 
# # # # 
# # # # 
# # # # 
# # # # # Below for backup purposes
# # # # ##############
# # # # # Final Filtering - Allele frequency and heterozygosity
# # # # ##############
# # # # 
# # # # 
# # # # 
# # # # # # Run post-filter diagnostic plots (really should make this step its own function or subscript or something)
# # # # # $TASSEL5 -SortGenotypeFilePlugin -inputFile $filtdir/1e_hemp_filter1.vcf.gz -outputFile $filtdir/1e_hemp_filter1.sorted.vcf.gz
# # # # # $TASSEL5 -vcf $filtdir/1e_hemp_filter1.sorted.vcf.gz -genotypeSummary site -export $filtdir/1e_sitesummary.filter1.txt
# # # # # $TASSEL5 -vcf $filtdir/1e_hemp_filter1.sorted.vcf.gz -genotypeSummary taxa -export $filtdir/1e_taxasummary.filter1.txt
# # # # # vcftools --gzvcf $filtdir/1e_hemp_filter1.sorted.vcf.gz  --site-depth --stdout > $filtdir/1e_site_depth.filter1.txt
# # # # # 
# # # # # # Plot 
# # # # # Rscript $scriptdir/1c_PlotGenoSummary.r --sitefile $filtdir/1e_sitesummary.filter1.txt --taxafile $filtdir/1e_taxasummary.filter1.txt --depthfile $filtdir/1e_site_depth.filter1.txt \
# # # # #    --restriction-frags $filtdir/1c_restriction_frag_lengths.txt --max-datapoints 10000 -o $filtdir/1e_hemp_filter1.summary_plots.png
# # # # # 
# # # # # 
# # # # #     
# # # # # # # Second round of filtering
# # # # # site_max_het=0.1
# # # # # site_min_maf=0.025
# # # # # Rscript $scriptdir/1d_GetFilterLists.r --sitefile $filtdir/1e_sitesummary.filter1.txt --taxafile $filtdir/1e_taxasummary.filter1.txt --depthfile $filtdir/1e_site_depth.filter1.txt \
# # # # #     --site-max-het $site_max_het --site-min-maf $site_min_maf \
# # # # #     --outtaxa $filtdir/1f_taxa_to_keep.filt2.txt --outsites $filtdir/1f_sites_to_keep.filt2.txt
# # # # # 
# # # # # # Perform actual filtering
# # # # # bcftools view --targets-file $filtdir/1f_sites_to_keep.filt2.txt --samples-file $filtdir/1f_taxa_to_keep.filt2.txt --output-type z --output-file $filtdir/1g_hemp_filter2.vcf.gz $filtdir/1e_hemp_filter1.vcf.gz
# # # # # 
# # # # # 
# # # # # # # Final post-filter diagnostic plots (mostly as a sanity check)
# # # # # $TASSEL5 -SortGenotypeFilePlugin -inputFile $filtdir/1g_hemp_filter2.vcf.gz -outputFile $filtdir/1g_hemp_filter2.sorted.vcf.gz
# # # # # $TASSEL5 -vcf $filtdir/1g_hemp_filter2.sorted.vcf.gz -genotypeSummary site -export $filtdir/1g_sitesummary.filter2.txt
# # # # # $TASSEL5 -vcf $filtdir/1g_hemp_filter2.sorted.vcf.gz -genotypeSummary taxa -export $filtdir/1g_taxasummary.filter2.txt
# # # # # vcftools --gzvcf $filtdir/1g_hemp_filter2.sorted.vcf.gz  --site-depth --stdout > $filtdir/1g_site_depth.filter2.txt
# # # # # 
# # # # # # Plot 
# # # # # Rscript $scriptdir/1c_PlotGenoSummary.r --sitefile $filtdir/1g_sitesummary.filter2.txt --taxafile $filtdir/1g_taxasummary.filter2.txt --depthfile $filtdir/1g_site_depth.filter2.txt \
# # # # #     --restriction-frags $filtdir/1c_restriction_frag_lengths.txt --max-datapoints 10000 -o $filtdir/1g_hemp_filter2.summary_plots.png
# # # # # 
# # # # # 
# # # # # # Copy filtered genotypes back to main directory
# # # # # cp $filtdir/1g_hemp_filter2.sorted.vcf.gz ./1g_hemp_filtered.sorted.vcf.gz
# # # # 
# # # # 
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
# # # # # ##############
# # # # # # Genomic prediction
# # # # # ##############
# # # # # 
# # # # # # Rename samples to match phenotype info (Tried with 'bcftools reheader', but it did the compression wrong somehow)
# # # # # source=1g_hemp_filtered.sorted.vcf.gz
# # # # # newfile=1h_hemp_filtered.sorted.renamed.vcf
# # # # # zcat $source | grep "^##" > $newfile   # Metadata lines
# # # # # zcat $source | grep "^#CHROM" | sed -r "s|1_bams_sorted/P.-HDGS-([0-9]+).aligned.sorted.bam|HDGS_\1|g" >> $newfile  # Header line
# # # # # zcat $source | grep -v "^#" >> $newfile   # Genotypes
# # # # # gzip $newfile
# # # # # 
# # # # # # Put phenotypes in TASSEL format; also do sqrt and log transformations, just in case
# # # # # Rscript $scriptdir/1i_FormatPhenosForTassel.r -i $hplc_results -o 1i_hdgs_phenos --log --sqrt
# # # # # 
# # # # # # Run genomic prediction with TASSEL
# # # # # $TASSEL5 -fork1 -vcf 1h_hemp_filtered.sorted.renamed.vcf.gz -ck -fork2 -t 1i_hdgs_phenos.tassel.txt -combine3 -input1 -input2 \
# # # # #     -GenomicSelectionPlugin -doCV true -kFolds 10 -nIter 1000 -endPlugin -export 1j_gs_results.txt
# # # # #     
# # # # # # Graph results as violin plots
# # # # # Rscript $scriptdir/1k_PlotPredictionAccuracy.r -i 1j_gs_results.txt -o 1k_gs_results
# # # # # 
# # # # # # Not great, unfortunately
