#! /usr/bin/Rscript

# Apply filters to determine which sites/taxa to keep

# Libraries
library(argparse)
library(RColorBrewer)
library(VennDiagram)

# Command-line arguments
options(stringsAsFactors=F)
parser=ArgumentParser()
parser$add_argument("-i", "--infiles", nargs="*", help="Input datasets (=first 2 columns of VCF file)")
parser$add_argument("-n", "--names", nargs="*", help="Names of input datasets")
parser$add_argument("-p", "--published", help="Sites from published dataset")
parser$add_argument("-o", "--outprefix", help="Output file prefix")
args=parser$parse_args()
# setwd('/home/jgwall/Projects/Misc/AddissuAyele_QuinoaGBS/2_CallSnps/')
# args=parser$parse_args(c("-i","2h_sites.filtered.txt", "2h_sites.highqual.txt", 
#                          "-n", "filt", "highqual",
#                          "-p", "2h_sites.published.txt",
#                          "-o", "99_tmp"))


# Load data
cat("Making Venn diagram of SNP overlaps\n")
pubdata = scan(args$published, what=character(), sep="\n")
indata = lapply(args$infiles, scan, what=character(), sep="\n")
names(indata) = args$names

# Force case the same
pubdata = toupper(pubdata)
indata = lapply(indata, toupper)

# Make colors
colors = brewer.pal(length(indata)+1, "Set2")
names(colors) = c("Published", names(indata))

# Plot
for(myname in names(indata)){
  mydata = list(pubdata, indata[[myname]])
  myfills = c(colors["Published"], colors[myname])
  outfile=paste(args$outprefix, myname, "png", sep=".")
  title=paste("Compare", myname,"to published set")
  venn.diagram(mydata, fill = myfills, alpha=0.5, cex=2, filename=outfile, 
               category.names=c("Published", myname), main=title,
               disable.logging=TRUE)
}