#! /usr/bin/Rscript

# Check and unify genotype and phenotype files
library(argparse)
library(tidyverse)
library(scales)
library(rrBLUP)

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

# Estimate heritabilities with rrBLUP
cat("\tEstimating heritabilities with rrBLUP\n")
phenos.rrblup = phenos %>%  # Reformat phenotypes
  rownames_to_column("taxa")
genos.rrblup = genos %>% # Reformat genotypes
  mutate(across(-taxa, ~rescale(.x, from=c(0,1), to=c(-1,1) ))) %>%
  column_to_rownames("taxa")
relations = A.mat(genos.rrblup) # Kinship matrix
models = lapply(names(phenos), function(mypheno){
  kin.blup(data=phenos.rrblup, geno="taxa", pheno=mypheno, K=relations)
})
herits = lapply(models, function(m){
  m$Vg / (m$Vg + m$Ve)
}) %>% unlist()
herits = data.frame(trait =names(phenos), h2_raw = herits,
                    h2=round(herits, digits=3))


# Plot phenotype distributions
phenodata = phenos %>%
  rownames_to_column("Taxa") %>%
  pivot_longer(cols=-Taxa, names_to="trait", values_to="value")
traitplot = ggplot(phenodata) +
  aes(x=value, fill=trait) +
  geom_histogram(bins=30) +
  facet_wrap(~trait, scales="free_x") +
  theme(legend.position="none") +
  labs(x="Trait value", y="Count", title="GWAS trait distributions") +
  geom_text(data=herits, mapping=aes(x=-Inf, y=Inf, label=paste("h2 =",h2)),
            hjust=0, vjust=1)
ggsave(traitplot, file=paste(args$outprefix, "_phenos.png", sep=""), width=8, height=8)


# Write out phenotypes
outphenos = data.frame(Taxa=rownames(phenos), phenos)
write.csv(outphenos, file=paste(args$outprefix, "_phenos.csv", sep=""), row.names = F)

# Convert genotypes to GAPIT range (0,1,2) and write out as RDS to save time later
# Confirm genos in 0-1 range
genorange = range(genos[,-1], na.rm=TRUE)
if(!identical(genorange, c(0,1))){
  stop("GENOTYPES NOT IN EXPECTED RANGE OF 0-1")
}
# Convert
genos.gapit = genos %>%
  mutate(across(-taxa, ~rescale(.x, from=c(0,1), to=c(0,2) )))
# Write out
write_rds(genos.gapit, file=paste(args$outprefix, "_genos.numeric.GAPIT.rds", sep=""))


