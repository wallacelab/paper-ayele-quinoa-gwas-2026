#! /usr/bin/Rscript

# This is a basic template to use for making R scripts
library(argparse)
library(tidyverse)

# Arguments
parser=ArgumentParser()
parser$add_argument("-a", "--prefile", help="MultiQC file pre filtering")
parser$add_argument("-b", "--postfile", help="MultiQC file post filtering")
parser$add_argument("-o", "--outprefix", help="Output file prefix")
args=parser$parse_args()
# setwd('/home/jgwall/Projects/Misc/AddissuAyele_QuinoaGBS/1_QualityControl/')  # Debug working directory
# args=parser$parse_args(c("-a","1a_FastQC.multiqc/multiqc_data/multiqc_general_stats.txt", 
#                          "-b","1b_Trimmed.multiqc/multiqc_data/multiqc_general_stats.txt",
#                          "-o",'99_tmp'))  # Debug test files

# Load data
pre=read.delim(args$prefile)
post=read.delim(args$postfile)

# Compile read retention stats
precount = pre %>%
  select(Sample, FastQC_total_sequences) %>%
  rename(pre_count=FastQC_total_sequences)
postcount = post %>%
  select(Sample, FastQC_total_sequences) %>%
  rename(post_count=FastQC_total_sequences) %>%
  mutate(Sample = sub(Sample, pattern="P$", repl="")) # Strip terminal "P"
counts = full_join(precount, postcount, by="Sample") %>%
  mutate(percent = post_count / pre_count)

# Write out stats
write.csv(counts, file=paste(args$outprefix, ".csv", sep=""), row.names=FALSE)

# Plot
plotdata = counts %>%
  pivot_longer(cols=-Sample, names_to="set", values_to="value") %>%
  mutate(set=factor(set, levels=c("pre_count", "post_count", "percent")))
myplots = ggplot(plotdata) +
  aes(x=Sample, y=value, fill=set) +
  geom_col() +
  facet_wrap(~set, nrow=1, scales="free_y") +
  theme(axis.text.x = element_blank(),
        legend.position="none") +
  labs(x="Samples", y="Read Count / Percent Retained")
ggsave(myplots, file=paste(args$outprefix, ".png", sep=""), width=8, height=3)
