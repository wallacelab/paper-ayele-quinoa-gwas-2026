#! /usr/bin/Rscript

# Check and unify genotype and phenotype files
library(argparse)
library(tidyverse)
library(openxlsx)

# Arguments
parser=ArgumentParser()
parser$add_argument("-g", "--infile", help="Input file of 2024 grenhouse data (csv)")
parser$add_argument("-p", "--physfile", help="Input file of 2025 physiological traits (xlsx)")
parser$add_argument("-s", "--seedfile", help="Input file of 2025 seed yield (xlsx)")
parser$add_argument("-u", "--updatefile", help="Input file of 2025 heights and corrected seed yields (xlsx)")
parser$add_argument("-o", "--outprefix", help="Output file prefix")
args=parser$parse_args()
# setwd('/home/jgwall/Projects/Misc/AddissuAyele_QuinoaGBS/0_Phenotypes')  # Debug working directory
# args=parser$parse_args(c("-g","2024_Greenhouse data.csv",
#                          "-p", "2025-Physiological traits.xlsx",
#                          "-s", "2025 Seed yield.xlsx",
#                          "-u", "2025_3TRT-Plant height _MTD_Yield2.xlsx",
#                          "-o",'99_tmp'))


# Load and process 2024 data
cat("Standardizing 2024 data\n")
raw2024=read.csv(args$infile) %>%
  mutate(Taxa = as.character(Taxa), ENTRY=as.character(ENTRY), TRT=as.character(TRT)) %>% # Make characters so treated as factors
  relocate(Taxa, TRT) %>%
  select(-ENTRY) %>%
  mutate(TRT = paste("TRT", TRT, sep="")) # So match 2025 and is a character vector


# Load & Process 2025 data - INCLUDES SPECIFICS OF WHICH SHEET TO LOAD
cat("Standardizing 2025 data\n")
raw.seed = read.xlsx(args$seedfile, sheet="Sheet1") %>%
  rename("Taxa" = `Plot#`, "Seedweight" = `Seed.wt..(g)`) %>%
  mutate(TRT=paste("TRT", TRT, sep="")) # Standardize with phys traits
# Read individual sheets (not combined)
phys.sheets = c("TRT1", "TRT2", "TRT3")
raw.phys.indiv = lapply(phys.sheets, function(s){
  mydata = read.xlsx(args$physfile, sheet=s)
  names(mydata) = sub(names(mydata), pattern="<Traits>", rep="<Trait>") # Fix inconsistency
  mydata %>% rename("Taxa" = "<Trait>") %>%
    mutate(TRT = s) %>%
    relocate(Taxa, TRT)
})
# Join together
metacols = c("Taxa", "TRT")
raw2025 = bind_rows(raw.phys.indiv) %>%
  left_join(raw.seed, by=metacols)


# Updated 2025 data with heights, maturity, and weights
cat("Standardizing updated 2025 data\n")
raw.update = lapply(phys.sheets, function(s){
  mydata = read.xlsx(args$updatefile, sheet=s)
  names(mydata) = sub(names(mydata), pattern="<Taxa>", rep="Taxa") # Fix taxa column
  names(mydata) = sub(names(mydata), pattern="Maturitydate", rep="Maturity.date") # Fix Maturity column inconsistency
  names(mydata) = sub(names(mydata), pattern="-.+", rep="") # Remove parts after hyphen
  names(mydata) = sub(names(mydata), pattern="\\.\\..+", rep="") # Remove parts after double period
  mydata %>% mutate(TRT = s) %>%
    rename(DTM=Maturity.date, Seedweight=Seed.wt) %>% # Standardize with other datasets
    relocate(Taxa, TRT)
}) %>% bind_rows()



# Write out
write.csv(raw2024, file=paste(args$outprefix, ".2024.csv", sep=""), row.names=FALSE)
write.csv(raw2025, file=paste(args$outprefix, ".2025.csv", sep=""), row.names=FALSE)
write.csv(raw.update, file=paste(args$outprefix, ".2025_update.csv", sep=""), row.names=FALSE)

