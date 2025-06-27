#! /usr/bin/Rscript

# Apply filters to determine which sites/taxa to keep

# Libraries
library(argparse)
library(RColorBrewer)
library(VennDiagram)

# Command-line arguments
options(stringsAsFactors=F)
parser=ArgumentParser()
parser$add_argument("-i", "--infiles", nargs=2, help="Input datasets (=first 2 columns of VCF file)")
parser$add_argument("-n", "--names", nargs=2, help="Names of input datasets")
parser$add_argument("-o", "--outprefix", help="Output file prefix")
args=parser$parse_args()
# setwd('/home/jgwall/Projects/Misc/AddissuAyele_QuinoaGBS/2_CallSnps/')
# args=parser$parse_args(c("-i","2h_sites.filtered.txt", "2h_sites.published.txt",
#                          "-n", "filt", "published",
#                          "-o", "99_tmp"))


# Load data
cat("Making Venn diagram of SNP overlaps\n")
indata = lapply(args$infiles, scan, what=character(), sep="\n")
names(indata) = args$names

# Force case the same
indata = lapply(indata, toupper)

# Make colors
colors = c("cornflowerblue", "firebrick")

# Plot
n1=names(indata)[1]
n2=names(indata)[2]
outfile=paste(args$outprefix, n1, n2, "png", sep=".")
title=paste("Compare", n1,"to",n2)
venn.diagram(indata, fill = colors, alpha=0.5, cex=2, filename=outfile, 
             category.names=names(indata), main=title,
             disable.logging=TRUE)
