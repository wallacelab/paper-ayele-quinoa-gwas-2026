#! /usr/bin/Rscript

# This is a basic template to use for making R scripts
library(argparse)
library(tidyverse)
library(ggpubr)

# Arguments
parser=ArgumentParser()
parser$add_argument("-i", "--infile", help="Input file of read join stats (sample, set, count)")
parser$add_argument("-o", "--outprefix", help="Output file prefix")
args=parser$parse_args()
# setwd('/home/jgwall/Projects/Misc/AddissuAyele_QuinoaGBS/1_QualityControl/')  # Debug working directory
# args=parser$parse_args(c("-i","1c_pairjoin_file_stats.txt",
#                          "-o",'99_tmp'))  # Debug test files

# Load data
counts = read.delim(args$infile)

# Pivot wider for some manipulations
wide = pivot_wider(counts, id_cols="sample", names_from="set", values_from="reads")

# Sanity check
mismatch = wide$forward != wide$reverse
if(any(mismatch, na.rm=TRUE)){
  problems = paste(wide$sample[mismatch], collapse=", ")
  warning("Forward and Reverse unassembled reads don't match for samples ", problems)
}

# Calculate percentage
wide$percent = wide$assembled / (wide$assembled + wide$forward + wide$discarded)

# Back to long format for plotting
plotdata = wide %>%
  select(-reverse) %>% # redundant
  rename(unassembled=forward, percent_joined=percent) %>%
  pivot_longer(cols=-sample, names_to="set", values_to="reads") %>%
  mutate(set = factor(set, levels=c("assembled", "unassembled", "discarded", "percent_joined")),
         sample = factor(sample)) 
  

# # Plot - Deprecated to split to keep axis scales together for top 2
# myplot = ggplot(plotdata) +
#   aes(x=sample, y=reads, fill=set) +
#   geom_col() +
#   facet_wrap(~set, nrow=2, scales="free_y") +
#   theme(axis.text.x = element_text(angle=90, vjust=0.5),
#         legend.position="none") +
#   ggtitle("Stats on read pair joining")
# ggsave(myplot, file=paste(args$outprefix, ".png", sep=""), width=5, height=5)

topdata = subset(plotdata, plotdata$set %in% c("assembled", "unassembled"))
bottomdata = subset(plotdata, plotdata$set %in% c("discarded", "percent_joined"))
topplot = ggplot(topdata) +
  aes(x=sample, y=reads, fill=set) +
  geom_col() +
  facet_wrap(~set, nrow=1) +
  theme(axis.text.x = element_blank(), legend.position="none") +
  scale_fill_manual(values=c("lightgreen", "coral")) +
  ggtitle("Read counts after joining")
bottomplot = ggplot(bottomdata) +
  aes(x=sample, y=reads, fill=set) +
  geom_col() +
  facet_wrap(~set, nrow=1, scales="free_y") +
  theme(axis.text.x = element_blank(), legend.position="none") +
  scale_fill_manual(values=c("darkred", "violet"))

myplot = ggarrange(topplot, bottomplot, nrow=2)
ggsave(myplot, file=paste(args$outprefix, ".png", sep=""), width=6, height=5)

# Write out text
write.csv(wide, file=paste(args$outprefix, ".csv", sep=""), row.names=FALSE)

