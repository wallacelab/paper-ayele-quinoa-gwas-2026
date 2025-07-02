#! /usr/bin/Rscript

# Test GAPIT for Addissu
library(GAPIT3)

setwd('/home/jgwall/Projects/Misc/AddissuAyele_QuinoaGBS/99_TestGAPIT/')

# Load data
phenos=read.csv("PhenoJ139.csv")
hmp = read.delim("../2_CallSnps/2f_snps.filtered.highqual.hmp.txt.gz", header=FALSE)
hmp.CD = read.delim("/home/jgwall/Projects/0_RawData/Sequencing/2024_12_13_Addisu_quinoa_GBS/CD_genomics_snp_calls_v2/Chenopodium_quinoa-dp4-miss0.2-maf0.05.hmp.gz", header=FALSE)

# Standardize names between them
phenos$Taxa = gsub(phenos$Taxa, pattern="2_CallSnps/2a_alignments/", repl="sample")
hmp[1,] = gsub(hmp[1,], pattern="2_CallSnps/2a_alignments/(.+).bam", repl="sample\\1")
hmp.CD[1,12:ncol(hmp.CD)] = paste("sample", as.character(hmp.CD[1,12:ncol(hmp.CD)]), sep="")

# Quick check; both of these should be equal to the number of samples (139 in this case)
sum(phenos$Taxa %in% hmp[1,])
sum(phenos$Taxa %in% hmp.CD[1,])


# Just do 1 trait to test
phenos = phenos[,1:2]


# Run GAPIT on Jason genos
dir.create("test_GAPIT_jason")
setwd("test_GAPIT_jason")
GAPIT(G=hmp, Y=phenos, model="GLM")
setwd("../")

# Run GAPIT on CD Genomics genos
dir.create("test_GAPIT_CDgenomics")
setwd("test_GAPIT_CDgenomics")
GAPIT(G=hmp.CD, Y=phenos, model="GLM")
setwd("../")