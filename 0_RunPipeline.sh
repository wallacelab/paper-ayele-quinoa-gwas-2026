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
pub_diversity=/home/jgwall/Projects/0_RawData/Diversity/Quinoa/Result/SNP/Cq_allpseudo_biallelic_minDP5_MaxMissing0.8_MAF0.01.Cq3Bsamples.variantsonly.vcf.gz
company_diversity=/home/jgwall/Projects/0_RawData/Sequencing/2024_12_13_Addisu_quinoa_GBS/CD_genomics_snp_calls_v2/Result/SNP/Chenopodium_quinoa-dp4-miss0.2-maf0.05.vcf.gz
srcdir=/home/jgwall/Projects/0_RawData/Sequencing/2024_12_13_Addisu_quinoa_GBS/01.rawdata/Chenopodium_quinoa-dp4-miss0.2-maf0.05.sorted.vcf.gz
adapterdir="/home/jgwall/Software/Python/miniconda3/envs/trimmomatic-0.39/share/trimmomatic-0.39-2/adapters/" # Location of trimmomatic adapters
rawdir=0_raw_data
num_cores=7  # Number of CPU cores to use


# Set up directories for major steps
scriptdir=0_Scripts
phenodir=0_Phenotypes
qcdir=1_QualityControl
trimdir=$qcdir/1b_Trimmed
joindir=$qcdir/1c_PairJoining
snpdir=2_CallSnps
gwasdir=3_GWAS
if [[ ! -e $rawdir ]]; then ln -s $srcdir $rawdir; fi  # Link raw data instead of copying
if [[ ! -e $qcdir ]] ; then mkdir $qcdir; fi
if [[ ! -e $snpdir ]] ; then mkdir $snpdir; fi
if [[ ! -e $gwasdir ]] ; then mkdir $gwasdir; fi

# Combine phenotypes for years 2024 and then 2025
# TODO - Compare 2025 physiology datasets to see what different/which to keep
# TODO - Figure out if anything needs to be combiend to blues
# TODO - Do QC on 2025 data
# TODO - Do programmatic cleaning of all data (outliers, boxcox transform, etc)
# TODO - Run GWAS on all data
bash $scriptdir/0_CleanPhenos.r $scriptdir $phenodir  # Note: a lot of needed files, etc. included in script instead of passing by argument


# # Basic quality control (fastqc, multiqc, adaptor trimming, etc. Looking at one of CD Genomics' BAM files, everything was still 144 bp, so no adapter trimming?)
# # NOTE: Not trimming off after a second cut site; will rely on chimera detection to filter those out
# adapters=$scriptdir/adapters.fa
# echo "" | cat $adapterdir/TruSeq3-PE-2.fa - $adapterdir/NexteraPE-PE.fa > $adapters  # Echo part is to put a newline between the files; using both to be conservative. Nextera transposon seqs definitely present.
# min_overlap=15  # Minimum overlap for joining paired ends reads with PEAR
# max_pval=0.001  # MAximum p-value for paired end joining with PEAR
# bash $scriptdir/1_QualityControl.sh $scriptdir $qcdir $rawdir $trimdir $pairdir $adapters $min_overlap $max_pval $num_cores


# # Call SNPs; supposedly mpileup only uses primary alignments, so that plus min quality 20 should remove chimeras
# max_depth=10000 # Maximum depth to pileup
# min_base_quality=20 # Minimum base quality for pileup
# subsample_for_stats=100000 # Subsample the raw SNP file to this amount to get genotyping stats
# minDP=5  # Minimum genotype depth to keep a genotype call
# maxDP=200 # Maximum genotype depth to keep a genotype call
# min_maf=0.05   # Minimum minor allele freq to keep a site
# max_missing=0.5 # Max fraction missing to keep a site
# max_het=0.10 # Max % heterozygous to keep a site
# max_mean_depth=55 # Max mean depth per site (to help remove pileups from paralogs)
# highqual_max_missing_site=0.2 # Max missing (sites) for high quality filters
# highqual_max_missing_taxon=0.25 # Max missing (taxon) for high quality filters
# highqual_min_mean_depth=20 # Minimum mean readdepth for high quality filters
# bash $scriptdir/2_CallSnps.sh $scriptdir $snpdir $trimdir $joindir $genome $max_depth $min_base_quality $subsample_for_stats \
#     $minDP $maxDP $min_maf $max_missing $max_het $max_mean_depth $highqual_max_missing_site $highqual_max_missing_taxon $highqual_min_mean_depth $pub_diversity \
#     $company_diversity $num_cores


# # Run GWAS with GAPIT 
# #genos=$scriptdir/mdp_genotype.short.hmp.txt.gz  # DEBUG Shorter version with 1k sites to debug quickly
# genos=$scriptdir/mdp_genotype.hmp.txt.gz  # Genotype file from Step 2, modified by coauthor to remove contigs and rename chromosomes
# phenos=$phenodir/mdp_traits.txt  # Phenotype file
# num_pcs=3   # Number of principal components to include
# num_perms=100   # Number of permutations to run to get null distribution of results
# num_resamples=100  # Number of resamples to run to get RMIPs
# resample_fraction=0.9  # Resample fraction number
# perm_p_cutoff="0.05 0.01"  # P-value cutoffs for permutations to count as null "hits"; make a single string with space delimits among different options
# rmip_cutoff="0.05" # RMIP cutoff for a hit to count as real; make a single string with space delimits among different options
# bash $scriptdir/3_RunGapit.sh $scriptdir $gwasdir $genos $phenos $num_pcs $num_perms $num_resamples $resample_fraction "$perm_p_cutoff" "$rmip_cutoff" $genome $num_cores


# TODO - Copy above and rerun with 2024 data. Will need to collapse to BLUPs first, but then rerun above just in a separate subfolder
