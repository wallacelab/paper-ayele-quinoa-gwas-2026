#! /usr/bin/Rscript

# Run permutation and resamplings to get RMIP for FARM-CPU
library(argparse)
library(GAPIT3)
library(dplyr)

# Arguments
parser=ArgumentParser()
parser$add_argument("-g", "--genofile", help="Input genotype file (RDS file, numeric matrix)")
parser$add_argument("-p", "--phenofile", help="Input file of phenotyes")
parser$add_argument("-t", "--trait", help="Which trait to run FamCPU on")
parser$add_argument("-c", "--covfile", help="Input file of covariates")
parser$add_argument("-k", "--kinshipfile", help="Input kinship matrix file")
parser$add_argument("--covariates", nargs="*", help="Which covariates to include")
parser$add_argument("-m", "--mapfile", help="Input GAPIT map file of genotypes")
parser$add_argument("-q", "--num-pcs", type="integer", default="5", help="How many principal components to include")
parser$add_argument("--permdir", help="Output directory for permutation results")
parser$add_argument("--resampledir", help="Output directory for resampling results")
parser$add_argument("--num-perms", type="integer", default="1", help="How many permutations to run of genotypes")
parser$add_argument("--num-resamples", type="integer", default=1, help="How many resamplings to run of genotypes")
parser$add_argument("-f", "--resample-fraction", type="double", default=0.8, help="What fraction of genotypes are included in each resampling")
args=parser$parse_args()
# setwd('/home/jgwall/Projects/Misc/AddissuAyele_QuinoaGBS/3_GWAS/')
# args=parser$parse_args(c("-g","3a_genos.numeric.rds",
#                          "-m","3a_genos.map.txt",
#                          "-p","3a_phenos.csv",
#                          "-k", "3a_kinship.txt",
#                          "-t", "GSW60",
#                          "-q", "3",
#                          "--permdir", "3b_Permutations/GSW60/",
#                          "--resampledir", "3d_Resampling/GSW60/",
#                          "--num-perms", "2",
#                          "--num-resamples", "2",
#                          "-f", "0.9"
# ))


cat("Running", args$num_perms,"FarmCPU permutations for", args$trait, "with", args$num_pcs, "PCs and covariates", args$covariates, "\n")
basedir=getwd()

# Load data
map = read.delim(args$mapfile)
genos=readRDS(args$genofile)
phenos=read.csv(args$phenofile)[,c("Taxa", args$trait)]
kinship=read.delim(args$kinshipfile, skip=3, row.names=1, header=FALSE)

# Format kinship for GAPIT
rownames(kinship) = basename(rownames(kinship)) # Fix taxa names
kinship = data.frame(taxa=rownames(kinship), kinship)

# Find shared taxa 
shared = intersect(genos$taxa, phenos$Taxa)

# Handle covariates if needed
if(!is.null(args$covfile)){
    covariates = read.csv(args$covfile)
    covariates = subset(covariates, select = names(covariates) %in% c("Taxa", args$covariates))
    
    # Remove null ones b/c missing covariate values cause obscure crash in GAPIT with PC calculation
    na_cov = is.na(covariates)
    covariates = subset(covariates, rowSums(na_cov) == 0)
    
    # Reduce to shared taxa
    shared = intersect(shared, covariates$Taxa)
    covariates = subset(covariates, covariates$Taxa %in% shared)
    
    # If no covariates left, set to null. Otherwise print out to verify
    if(ncol(covariates)==1){
        covariates = NULL
    }else{
        cat("\n## Including covariates",names(covariates)[-1], "##\n")
    }
}else{
    covariates = NULL
}

# Cut genos & phenos down to shared taxa.
genos = subset(genos, genos$taxa %in% shared)
phenos = subset(phenos, phenos$Taxa %in% shared)

# Impute missing data to mean for each site (otherwise GAPIT crashes)
#genos[is.na(genos)] = 0
missing_to_mean = function(x){
  x[is.na(x)] = mean(x, na.rm=TRUE)
  return(x)
}
genos = genos %>%
  mutate(across(-taxa, missing_to_mean))



# # Generate kinship matrix (stays constant) - DEPRECATED
# kinship = GAPIT3:::GAPIT.kinship.VanRaden(genos)
# kinship = data.frame(taxa=rownames(kinship), kinship, check.names=FALSE)

# Calculate number of genotypes per resample
genos_per_resample = args$resample_fraction * nrow(genos)

# Run permutations
for(i in 1:args$num_perms){
    cat("\n#### Running permutation",i,"####\n\n")
    set.seed(i) # Set random seed for reproducibility
    
    # Set up working directory
    outdir=paste(args$permdir, "/perm", i, sep="")
    dir.create(outdir, recursive=TRUE)
    setwd(outdir)
    
    # Scramble genotypes among taxa
    permutations = sample(1:nrow(genos), replace=F) 
    genos.perm = data.frame(taxa=genos$taxa, genos[permutations,-1])
    
    # Subset genotypes so have as many as will be in each resample
    to_keep = sample(1:nrow(genos.perm), replace=F, size=genos_per_resample)
    genos.perm = genos.perm[to_keep,]
    
    # # Run FarmCPU
    results = GAPIT(Y=phenos, GD=genos.perm, GM=map, CV=covariates, KI=kinship,
          PCA.total=args$num_pcs, model=c("FarmCPU"), file.output=T)

    # Reset to original directory
    setwd(basedir)
}



# Run resamples
for(i in 1:args$num_resamples){
    cat("\n#### Running resample",i,"for",args$trait,"####\n\n")
    set.seed(i) # Set random seed for reproducibility
    
    # Set up working directory
    outdir=paste(args$resampledir, "/resample", i, sep="")
    dir.create(outdir, recursive=TRUE)
    setwd(outdir)
    
    # Subset genotypes
    to_keep = sample(1:nrow(genos), replace=F, size=genos_per_resample)
    subgenos = genos[to_keep,]
    
    # Run FarmCPU
    results = GAPIT(Y=phenos, GD=subgenos, GM=map, CV=covariates, KI=kinship,
                    PCA.total=args$num_pcs, model=c("FarmCPU"), file.output=T)
    
    # Reset to original directory
    setwd(basedir)
}