#! /usr/bin/Rscript

# Plot GAPIT results as Manhattan-like plots

library(argparse)
library(tidyverse)
library(ggthemes)
library(ggpubr)
library(scales)


# Arguments
parser=ArgumentParser()
parser$add_argument("-i", "--infiles", nargs="*", help="Input files of parsed RMIP values from step 3e_")
parser$add_argument("-w", "--winsize", default=5e6, type="double", help="Window size (in bp) for RMIP sliding window")
parser$add_argument("-s", "--step", default=5e5, type="double", help="Step size (in bp) for RMIP sliding window")
parser$add_argument("--drop-contigs", default=FALSE, action="store_true", help="Whether to drop contigs from the plot")
parser$add_argument("--index", help="Samtools index of chromosome sizes")
parser$add_argument("--rmip-cutoff", type='double', default=0.05, help="RMIP cutoff to be consider a genuine hit")
parser$add_argument("-o", "--outprefix", help="Output file prefix")
args=parser$parse_args()
# setwd('/home/jgwall/Projects/Misc/AddissuAyele_QuinoaGBS/3_GWAS/')
# args=parser$parse_args(c("-i","3e_Results.FarmCPU/3e_Seedweight.p0.01.rmip0.05.csv",
#                        "3e_Results.FarmCPU/3e_GSW56.p0.01.rmip0.05.csv",
#                        "3e_Results.FarmCPU/3e_GSW60.p0.01.rmip0.05.csv",
#                        "--index", "3f_genome_index.txt",
#                        "--drop-contigs",
#                        "-o", "99_tmp"))

# Load data
index=read.delim(args$index, header=FALSE,
                 col.names=c("chr", "length","x1", "x2", "x3"))
rmips = lapply(args$infiles, function(infile){
  farm = read.csv(infile)
})

# Assign trait names to rmip data
traits = basename(args$infiles) %>%
  sub(pattern="3e_", repl="") %>%
  sub(pattern=".p.+rmip.+.csv", rep="")
names(rmips) = traits

# Convert index into a set of chromosome offsets for plotting
if(args$drop_contigs){
  index = subset(index, !grepl("contig", index$chr, ignore.case = TRUE))
}
offsets = index %>%
  select(chr, length) %>%
  arrange(chr) %>%
  mutate(offset=cumsum(length)) %>%
  mutate(offset=c(0, offset[1:(n()-1)])) %>% # Make sure starts at 0
  mutate(offset = offset + 1:n() * 5e6) %>% # Add padding between chromosomes
  mutate(midpoint = offset + length/2) # Calculate midpoint of each chromosome)

# Helper function to get x position for plotting a point given a vector of chroms and positions
get_offset = function(chr, pos){
  myoffset = offsets[as.character(chr),"offset"]
  return(myoffset + pos)
}


# Plot a single manhattan plot
plot_manhattan = function(trait_name){
  
  # Load in data
  cat("Plotting data for",trait_name, "\n")
  rmip_data = rmips[[trait_name]]
  
  # If no data, return null plot
  if(nrow(rmip_data)==0){
    nullplot = ggplot() + 
      annotate("text", x = 4, y = 25, size=8,
               label = paste("No hits for",trait_name), color="red") + 
      theme_void()  
    return(nullplot)
  }
  colorkey = c(yes="red", no="black")
  
  # Plot raw RMIPs
  myplot = ggplot(rmip_data) +
    aes(x = get_offset(chr, pos), y=rmip, color=significant, fill=significant) + 
    geom_point(size=2, shape=25) +
    scale_fill_manual(values=colorkey) +
    scale_color_manual(values=colorkey)+ 
    ylim(c(0,1)) +
    theme_few() + # Low-clutter theme
    scale_x_continuous(breaks=offsets$midpoint, minor_breaks = offsets$offset, labels=rownames(offsets)) +  # Set labels and line positions for X axis
    theme(legend.position = "none") +  # Turn off legend
    labs(x="Chromosome", y="RMIP score") + # X and Y axis labels
    ggtitle(trait_name) +
    theme(title = element_text(size=10, face='bold'))
  

  # Cumulative Farm-CPU distribution
  cum_color='gray35'
  for(mychrom in rownames(offsets)){
    
    # Get cumulative RMIP values
    #cat("\tOn chrom",mychrom,"with length",offsets[mychrom, "length"], "\n") # For debuggin
    mypoints = data.frame(chr=mychrom, pos=seq(0, offsets[mychrom, "length"], by=args$step), cum_rmip=0,
                          significant=NA) # "significant" column s ggplot doesn't throw an error
    for(i in 1:nrow(mypoints)){
      mypos = mypoints[i, "pos"]
      nearby_hits = subset(rmip_data, rmip_data$chr==mychrom & 
                             abs(rmip_data$pos - mypos) <= args$winsize/2)
      mypoints$cum_rmip[i] = sum(nearby_hits$rmip)
    }
    
    # Add to plot
    myplot = myplot +
      geom_line(data=mypoints, mapping=aes(y=cum_rmip), color=cum_color,  linewidth=0.35) 
  }
  
  # Add horizontal reference line
  myplot = myplot +
    geom_hline(yintercept=args$rmip_cutoff, color="red", linewidth=0.35, linetype="dashed")

  return(myplot)
}

manhattans = lapply(traits, plot_manhattan)


# Arrange all plots and save
outplot = ggarrange(plotlist=manhattans, ncol=1)
ggsave(outplot, file=paste(args$outprefix, ".png", sep=""),
       width=8, height=2*length(traits))

