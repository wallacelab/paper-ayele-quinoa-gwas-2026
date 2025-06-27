#! /usr/bin/Rscript

# This is a basic template to use for making R scripts
library(argparse)
library(tidyverse)

# Arguments
parser=ArgumentParser()
parser$add_argument("-g", "--genofile", help="Input genotype file (numeric format from TASSEL)")
parser$add_argument("-p", "--phenofile", help="Input files of phenotyes")
parser$add_argument("-o", "--outprefix", help="Output file prefix")
parser$add_argument("-c", "--covariates", nargs="*", help="Which traits are the covariates")
args=parser$parse_args()
# setwd('/home/jgwall/Projects/Misc/AddissuAyele_QuinoaGBS/')  # Debug working directory
# args=parser$parse_args(c("-p","0_Scripts/mdp_traits.txt",
#                          "-g",'3_GWAS/3a_genos.numeric.txt',
#                          "-o",'99_tmp'))


# Load data
cat("Loading and comparing genetic and phenotype data\n")
phenos = read.delim(args$phenofile, row.names=1)
genos = read.delim(args$genofile, skip=1, row.names=1,
                   colClasses="character", check.names=FALSE) # Load as character for speed

# Convert genos to numeric
genos = genos %>%
  mutate(across(everything(), as.numeric))

# Fix genotype names by removing path information
genos = genos %>% rownames_to_column("taxa")
genos$taxa = basename(genos$taxa)


# Compare genos and phenos (so make sure have taxa the same)
display = function(a, b, text){
    percent = round(length(a) / length(b) * 100, digits=1)
    cat("\t",length(a),"out of",length(b),text, "taxa are shared (",percent,"% )\n" )    
}
shared = intersect(rownames(phenos), genos$taxa)
display(shared, rownames(phenos), "trait")
display(shared, genos$taxa, "genotype")


# Plot phenotype distributions
phenodata = phenos %>%
  rownames_to_column("Taxa") %>%
  pivot_longer(cols=-Taxa, names_to="trait", values_to="value")
traitplot = ggplot(phenodata) +
  aes(x=value, fill=trait) +
  geom_histogram(bins=30) +
  facet_wrap(~trait, scales="free_x") +
  theme(legend.position="none") +
  labs(x="Trait value", y="Count", title="GWAS trait distributions")
ggsave(traitplot, file=paste(args$outprefix, "_phenos.png", sep=""), width=8, height=8)

# Write out phenotypes
outphenos = data.frame(Taxa=rownames(phenos), phenos)
write.csv(outphenos, file=paste(args$outprefix, "_phenos.csv", sep=""), row.names = F)

# Write out genotypes as RDS to save time later
write_rds(genos, file=paste(args$outprefix, "_genos.numeric.rds", sep=""))
