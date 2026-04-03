#! /usr/bin/Rscript

# Reformat cleaned phenotypes for GAPIT

# Load libraries
library(argparse)
library(tidyverse)

# Arguments
parser=ArgumentParser()
parser$add_argument("-i", "--infile", help="Input file of standardized pheno data (csv)")
parser$add_argument("-p", "--priority-traits", nargs="*", help="Strings indicating the start of priority phenotypes (used for splitting outfiles)")
parser$add_argument("-o", "--outprefix", help="Output file prefix")
args=parser$parse_args()
# setwd('/home/jgwall/Projects/Misc/AddissuAyele_QuinoaGBS/0_Phenotypes/')  # Debug working directory
# args=parser$parse_args(c("-i","0d_phenos.2025_update.clean.filtered.csv",
#                          "-p",'GSW', 'PhiPS',
#                          "-o",'99_tmp'))

cat("Formatting phenotypes in", args$infile,"for GAPIT\n")
phenos = read.csv(args$infile)

# Get list of traits
metacols = c("Taxa", "TRT") # Metadata columns
traits = names(phenos) %>% 
  setdiff(metacols) %>%
  sort()

# Identify target traits to port over
targets = lapply(args$priority_traits, function(mytrait){
  pattern = paste("^", mytrait, sep="") # Make sure matches as beginning
  matches = grep(traits, pattern=mytrait)
  return(traits[matches])
})

# Print out diagnostics to check
for(i in 1:length(targets)){
  cat("\tTarget",args$priority_traits[i],"matches",targets[[i]],"\n")
}
priorities = unlist(targets) # Make into a single vector
notpicked = setdiff(traits, priorities)
cat("\tTraits not prioritized:",notpicked,"\n")


# Split by priorities
phenos.top = phenos[,c(metacols, priorities)]
phenos.bottom = phenos[,c(metacols, notpicked)]

# Combine trait and treatment
trait_plus_treatment = function(mydata){
  mytraits = setdiff(names(mydata), metacols)
  if(length(mytraits)==0){
    return(data.frame(Taxa=NA))
  }
  mydata = split(mydata, mydata$TRT)
  mydata = lapply(mydata, function(m){
    trt = unique(m$TRT)
    m = rename_with(m, .cols = all_of(mytraits), .fn = ~ paste(.x, trt, sep="_")) %>%
      select(-TRT)
  }) %>%
    reduce(left_join, by="Taxa")
  return(mydata)
}
phenos.top = trait_plus_treatment(phenos.top)
phenos.bottom = trait_plus_treatment(phenos.bottom)


# Format for GAPIT and write out
output_gapit = function(mydata, outfile){
  mydata = mydata %>% rename("<Taxa>" = "Taxa")
  write.table(mydata, file=outfile, row.names=FALSE,
              col.names=TRUE, sep='\t', quote=FALSE)
}
output_gapit(phenos.top, outfile=paste(args$outprefix, ".priorities.gapit.txt", sep=""))
output_gapit(phenos.bottom, outfile=paste(args$outprefix, ".remainder.gapit.txt", sep=""))




