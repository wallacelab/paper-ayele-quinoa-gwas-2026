#! /usr/bin/Rscript

# Make Manhattan plots using GLM and RMIP data

library(argparse)
library(tidyverse)
library(ggthemes)
library(scales)
library(gridExtra)
library(ggpubr)
library(grid)
library(colorspace)
#library(patchwork)
library(png)
library(gridExtra)
library(grid)

# Arguments
parser=ArgumentParser()
# RMIP arguments
parser$add_argument("-r", "--rmip-infiles", nargs="*", help="Input files of parsed RMIP values from step 3e_")
parser$add_argument("-w", "--winsize", default=5e6, type="double", help="Window size (in bp) for RMIP sliding window")
parser$add_argument("-s", "--step", default=5e5, type="double", help="Step size (in bp) for RMIP sliding window")
parser$add_argument("--rmip-cutoff", type='double', default=0.05, help="RMIP cutoff to be consider a genuine hit")
# GLM arguments
parser$add_argument("-g", "--glm-infiles", nargs="*", help="Input files of GLM results values from step 3h_")
parser$add_argument("-p", "--pval-cutoff", type="double", default=0.01, help="Cutoff for empirical p-value to be considered significant")
# Other arguments
# parser$add_argument("--index", help="Samtools index of chromosome sizes")
parser$add_argument("-l", "--legend-file", help="File with custom legend to add below multi-treatment composite plots")
parser$add_argument("--drop-contigs", default=FALSE, action="store_true", help="Whether to drop contigs from the plot")
parser$add_argument("--debug", default=FALSE, action="store_true", help="Debugging (to make some steps faster")
parser$add_argument("-o", "--outprefix", help="Output file prefix")
args=parser$parse_args()
# setwd('/home/jgwall/Projects/Misc/AddissuAyele_QuinoaGBS/3_GWAS/Priorities2025/')
# args=parser$parse_args(c("-r","3e_Results.BLINK/3e_PhiPS60_TRT1.p0.05.rmip0.05.csv",
#                        "3e_Results.BLINK/3e_PhiPS60_TRT2.p0.05.rmip0.05.csv",
#                        "-g","3g_TASSEL_GLM/3h_glm_results_10.sites.txt.gz",
#                        "3g_TASSEL_GLM/3h_glm_results_28.sites.txt.gz",
#                        "-l", "../../0_Scripts/manhattan-legend.png",
#                        "--debug",
#                        "-o", "99_tmp"))


# # Output figure size (in inches)
# out.width=8
# #out.height=3 # Depends on the specific plot
# out.dpi=600
# legend.file="wallace-legend.png"

colorscale=c("#D81B60", "#1E88E5", "#FFC107")

##########
# Load data
##########

cat("Combining manhattan plots for", length(args$rmip_infiles), "RMIP files and",
    length(args$glm_infiles), "GLM files\n")

# Load RMIP data
load_rmip = function(infile){
  dummy_data = 
  if(file.exists(infile)){
    cat("\tLoading RMIP results from",infile,"\n")
    rmip = read.csv(infile)
  }else{
    cat("\tUnable to locate RMIP file",infile,"; creating dummy data\n")
    rmip = data.frame(snp=character(), chr=integer(), pos=integer(), rmip=numeric(), significant=character())
  }
  names(rmip)[names(rmip)=="chr"] = "Chr"
  names(rmip)[names(rmip)=="pos"] = "Pos"

  return(rmip)
}
rmips = lapply(args$rmip_infiles, load_rmip)
names(rmips) = basename(args$rmip_infiles) %>% # Assign trait names
  sub(pattern="3e_", repl="") %>%
  sub(pattern=".p.+rmip.+.csv", rep="")

# Load GLM data
load_glm = function(infile, debug=FALSE){
  if(!is.null(infile) && file.exists(infile)){
    cat("\tLoading GLM results from",infile,"\n")
    glm = read.delim(infile)
    
    # For debugging, take random sample of GLM results to speed process
    if(debug){
      set.seed(1);
      glm = glm[sample(1:nrow(glm), size=1000, replace=F),]
    }
  }else {
    cat("\tUnable to locate GLM file",infile,"; creating dummy data\n")
    glm = data.frame(Chr=1, Pos=1, p=c(1, 0.1), perm_p=1)
  }
  glm$log.p = -log10(glm$p)  # -log10 of p-value
  return(glm)
}
glms = lapply(args$glm_infiles, load_glm, debug=args$debug)
names(glms) = sapply(glms, function(x){unique(x$Trait)})

############
# Chromosome offsets
############

# Helper function to get max position per chromosome
get_max = function(x){
  if(nrow(x)==0){  # To handle empty RMIP results
    return(data.frame(Chr=1, Pos=1))
  }
  x %>%
    select(Chr, Pos) %>%
    group_by(Chr) %>%
    summarize(Pos=max(Pos))
}
# Make an index of chromosome lengths based on max positions
glm_max = lapply(glms, get_max) %>% bind_rows()
rmip_max = lapply(rmips, get_max) %>% bind_rows()
index =  bind_rows(glm_max, rmip_max) %>%
  get_max() %>%
  rename(length=Pos)

# Calculate chromosome offsets
offsets = data.frame(row.names = index$Chr, length=index$length)
offsets$offset = c(0, cumsum(offsets$length)[1:(nrow(offsets)-1)] ) +
  5e6 * 1:nrow(offsets)  # 5e6 = Padding between chromosomes
offsets$midpoint = offsets$offset + (offsets$length/2) # Calculate midpoint of each chromosome

# Helper function to get x position for plotting a point given a vector of chrom and positions
get_offset = function(chr, pos){
  myoffset = offsets[as.character(chr),"offset"]
  return(myoffset + pos)
}

############
# Plot Manhattans (1 per treatment)
############

# Plot a single manhattan plot
plot_manhattan = function(trait_name){

    # Load in data
    cat("\tPlotting data for",trait_name, "\n")
    glm = glms[[trait_name]]
    rmip = rmips[[trait_name]]
    
    if(is.null(glm)){
      cat("Cannot plot due to null results for", trait_name,"GLM\n")
      return(NULL)
    }
    if(is.null(rmip)){
      cat("Cannot plot due to null results for", trait_name,"RMIP\n")
      return(NULL)
    }
    
    # Subset out permutation significance
    glm.sig = subset(glm, glm$perm_p <= args$pval_cutoff)
    
    # Reformat RMIP data
    rmip.all = rmip # Backup of complete RMIP data
    rmip = rmip %>% filter(rmip >= args$rmip_cutoff)
    
    # Set up plotting colors
    chrom_colors=RColorBrewer::brewer.pal(3, "Set1")[1:2]
    my_alpha = 0.5
    if(all(glm$Pos==1)){ # Set alpha to 0 for GLM if is a dummy dataset (all positions at 1)
        my_alpha=0
    }
 
    # Basic plotting
    glm_offset = floor(min(glm$log.p))  # Offset by nearest integer below
    myplot = ggplot(glm) +
        aes(x=get_offset(Chr, Pos), y=log.p - glm_offset, color=(Chr%%2 == 0)) + 
        geom_point(alpha=my_alpha) +  # Make a scatterplot
        scale_color_manual(values=chrom_colors) +
        theme_few() + # Low-clutter theme
        scale_x_continuous(breaks=offsets$midpoint, minor_breaks = offsets$offset, labels=rownames(offsets)) +  # Set labels and line positions for X axis
        theme(legend.position = "none") +  # Turn off legend
        #labs(x="Chromosome", y=expression(bold(paste("-", log[10], " p-value (GLM)")))) + # X and Y axis labels
        labs(x=element_blank(), y=element_blank()) + # X and Y axis labels
        ggtitle(trait_name) +
        theme(title = element_text(size=10, face='bold')) # Set title 
    
    # Add significant GLM hits
    myplot = myplot +
        geom_point(data=glm.sig, size=2, color='black', fill=alpha('white', 0), shape=21, stroke=1.5)
    
    # Add right Y axis to plot and adjust left axis to reflect GLM offset
    if(nrow(rmip.all)>0){
        rmip.max = max(c(args$rmip_cutoff, rmip.all$rmip))
    }else{
        rmip.max=1
    }
    scale_trans = rmip.max / max(glm$log.p)    
    myplot = myplot +
        scale_y_continuous(labels = function(x){x+glm_offset},
          sec.axis = sec_axis(transform = ~. * scale_trans, name=element_blank()))
    
    # Cumulative rmip-CPU distribution
    if(nrow(rmip.all) > 0){
        cum_color='gray35'
        for(mychrom in rownames(offsets)){
            
            # Get cumulative RMIP values
            mypoints = data.frame(Chr=mychrom, Pos=seq(0, offsets[mychrom, "length"], by=args$step), cum_rmip=0)
            for(i in 1:nrow(mypoints)){
                mypos = mypoints[i, "Pos"]
                nearby_hits = subset(rmip.all, rmip.all$Chr==mychrom & 
                                         abs(rmip.all$Pos - mypos) <= args$winsize/2)
                mypoints$cum_rmip[i] = sum(nearby_hits$rmip)
            }
            
            # Get translated Y values
            mypoints$rmip.plot = rescale(mypoints$cum_rmip, from=c(0, max(rmip.all$rmip)), to=c(0, max(glm$log.p)))
            
            # Add to plot
            myplot = myplot +
                geom_line(data=mypoints, mapping=aes(y=rmip.plot), color=cum_color, linewidth=0.35) 
        }
        
        # Add horizontal reference line
        ref_line_y = rescale(args$rmip_cutoff, from=c(0, rmip.max), to=c(0, max(glm$log.p)))
        myplot = myplot +
            geom_hline(yintercept=ref_line_y, color=cum_color, size=0.35, linetype="dashed")
    }
    
    # Significant rmipCPU hits by RMIP
    if(nrow(rmip) > 0){
        rmip$rmip.plot = rescale(rmip$rmip, from=c(0, max(rmip$rmip)), to=c(0, max(glm$log.p)))
        rmip_color="black"
        myplot = myplot +
            #geom_segment(data=rmip, mapping=aes(xend=get_offset(Chr, Pos), y=0, yend=rmip.plot), color=rmip_color, size=0.5) +
            geom_point(data=rmip, mapping=aes(y=rmip.plot), size=2, color=rmip_color, fill=rmip_color, shape=25) 
    }
    
    return(myplot)
}


# Helper function to make plots for each phenotype
get_plots = function(traits){
  # Make plots
  myplots=list()
  for(i in 1:length(traits)){
    myplots[[i]] = plot_manhattan(traits[i])
  }

  # Arrange into grid
  chrom_label=text_grob("Chromosome", size = 12, face = "bold")
  left_label =text_grob(expression(bold(paste("-", log[10], " p-value (GLM)"))), size = 12, rot=90)
  right_label =text_grob("RMIP", size = 12, face='bold', rot=-90)
  #legend=rasterGrob(readPNG(legend.file))
  grid.arrange(arrangeGrob(grobs=myplots, ncol=1), #legend,
               nrow=2, heights=c(10,1), bottom=chrom_label, left=left_label, right=right_label)
}


##################
# Plot Manhattans - Combined across treatmetns (color differentiating)
##################

# Note: Lots of reused code here. Prefer to keep it separate so can tweak the two formats separately

# Make a manhattan plot combined across years
get_combined_plot = function(traits){
  
  # Load in data
  cat("\tPlotting single manhattan for",traits, "\n")
  trait_name = traits %>%
    sub(pattern="_TRT.+", repl="")  %>%
    unique()
  
  # Format and combined data
  glm = lapply(traits, function(t){
    glms[[t]] %>%
      mutate(treatment=sub(Trait, pattern=".+_(TRT.)", repl="\\1"))
  }) %>%
    bind_rows()
  rmip = lapply(traits, function(t){
    rmips[[t]] %>%
      mutate(Trait=t,
        treatment=sub(Trait, pattern=".+_(TRT.)", repl="\\1"),
        snp=as.character(snp))  # to fix subtle bug
  }) %>%
    bind_rows()
  
  # Sanity checks
  if(is.null(glm)){
    cat("Cannot plot due to null results for", trait_name,"GLM\n")
    return(NULL)
  }
  if(is.null(rmip)){
    cat("Cannot plot due to null results for", trait_name,"RMIP\n")
    return(NULL)
  }
  
  # Subset out permutation significance
  glm.sig = subset(glm, glm$perm_p <= args$pval_cutoff)
  
  # Reformat RMIP data
  rmip.all = rmip # Backup of complete RMIP data
  rmip = rmip %>% filter(rmip >= args$rmip_cutoff)
  
  # Set up plotting colors
  treat_levels = unique(c(glm$treatment, rmip$treatment)) %>% sort()
  #chrom_colors=RColorBrewer::brewer.pal(3, "Set1")[1:length(treat_levels)]
  chrom_colors=colorscale[1:length(treat_levels)] 
  names(chrom_colors) = treat_levels
  #chrom_colors_dark = chrom_colors %>% darken(amount=0.5)
  # my_alpha = 0.5
  # if(all(glm$Pos==1)){ # Set alpha to 0 for GLM if is a dummy dataset (all positions at 1)
  #   my_alpha=0
  # }
  
  # Basic plotting
  glm_offset = floor(min(glm$log.p))  # Offset by nearest integer below
  myplot = ggplot(glm) +
    aes(x=get_offset(Chr, Pos), y=log.p - glm_offset, color=treatment) + 
    #geom_point(alpha=0.05) + 
    geom_point(mapping=aes(color=treatment, color = after_scale(darken(color, 0.3))), alpha=0.05) +  # Make a scatterplot with darker colors
    scale_color_manual(values=chrom_colors) +
    theme_few() + # Low-clutter theme
    scale_x_continuous(breaks=offsets$midpoint, minor_breaks = offsets$offset, labels=rownames(offsets)) +  # Set labels and line positions for X axis
    labs(x="Position", y=expression(bold(paste("-", log[10], " p-value (GLM)"))),
         color="Treatment") + # X and Y axis labels
    ggtitle(trait_name) +
    theme(title = element_text(size=10, face='bold')) + # Set title 
    guides(alpha="none",  # Remove alpha transparency legend
           color = guide_legend(override.aes = list(alpha = 1, size=4)))  # Make legend points non-transparent

  # Add significant GLM hits
  myplot = myplot +
    geom_point(data=glm.sig, size=2, color="white", shape=21, stroke=2, show.legend=FALSE) + # White stroke
    geom_point(data=glm.sig, size=2, fill=alpha('white', 0), shape=21, stroke=1.5, show.legend=FALSE)
  
  # Add right Y axis to plot and adjust left axis to reflect GLM offset
  if(nrow(rmip.all)>0){
    rmip.max = max(c(args$rmip_cutoff, rmip.all$rmip))
  }else{
    rmip.max=1
  }
  scale_trans = rmip.max / max(glm$log.p)    
  myplot = myplot +
    scale_y_continuous(labels = function(x){x+glm_offset},
                       sec.axis = sec_axis(transform = ~. * scale_trans, name="RMIP"))

  # Add RMIPs (if exist)
  if(nrow(rmip.all) > 0){
    for(mytreat in treat_levels){
      cum_color=chrom_colors[mytreat]
      for(mychrom in rownames(offsets)){

        # Get cumulative RMIP values
        mypoints = data.frame(Chr=as.numeric(mychrom), Pos=seq(0, offsets[mychrom, "length"], by=args$step), cum_rmip=0)
        for(i in 1:nrow(mypoints)){
          mypos = mypoints[i, "Pos"]
          nearby_hits = subset(rmip.all, rmip.all$Chr==mychrom & rmip.all$treatment==mytreat &
                                 abs(rmip.all$Pos - mypos) <= args$winsize/2)
          mypoints$cum_rmip[i] = sum(nearby_hits$rmip)
        }

        # Get translated Y values
        mypoints$rmip.plot = rescale(mypoints$cum_rmip, from=c(0, max(rmip.all$rmip)), to=c(0, max(glm$log.p)))

        # Add cumulative RMIP trace to plot
        myplot = myplot +
          geom_line(data=mypoints, mapping=aes(y=rmip.plot), color="white", linewidth=1, alpha=0.5) +
          geom_line(data=mypoints, mapping=aes(y=rmip.plot), color=cum_color, linewidth=0.5)
      }

    }
    # Add horizontal reference line
    ref_line_y = rescale(args$rmip_cutoff, from=c(0, rmip.max), to=c(0, max(glm$log.p)))
    myplot = myplot +
      geom_hline(yintercept=ref_line_y, color="black", linewidth=0.35, linetype="dashed")
  }
  
  # Significant rmipCPU hits by RMIP
  if(nrow(rmip) > 0){
    rmip$rmip.plot = rescale(rmip$rmip, from=c(0, max(rmip$rmip)), to=c(0, max(glm$log.p)))
    #rmip_color="black"
    myplot = myplot +
      geom_point(data=rmip, mapping=aes(y=rmip.plot, fill=treatment), color="black", stroke=0.5, size=2.5, shape=25, show.legend = FALSE) +
      scale_fill_manual(values=chrom_colors)
  }
  
  # Add legend
  if(!is.null(args$legend_file)){
    legend = readPNG(args$legend_file, native=TRUE) %>%
      rasterGrob
    myplot = grid.arrange(arrangeGrob(grobs=list(myplot), ncol=1), legend, nrow=2, heights=c(10,1))
    
  }
  
  return(myplot)
  
}


###############
# Actual plotting section
###############



# Plot trait groupings by trait
stems = c(names(glms), names(rmips)) %>%
  sub(pattern="_TRT[0-9]", repl="") %>%  # Group 3 treatments
  unique() %>%
  sort()
# Loop over stems
for(stem in stems){
  
  cat("Making plots based on trait stem",stem,"\n")
  traits = union(names(glms), names(rmips))
  traits = traits[grep(traits, pattern=paste("^", stem, sep=""))] %>%
    sort()
  
  # Plot and save (separate for each treatment)
  myplots = get_plots(traits)
  outfile1=paste(args$outprefix, stem, "png", sep=".")
  ggsave(myplots, file=outfile1, width=12, height=2*length(traits))
  
  # Plot and save (combined across TRTs in one image)
  combined_plot = get_combined_plot(traits)
  outfile2=paste(args$outprefix, stem, "combined", "png", sep=".")
  ggsave(combined_plot, file=outfile2, width=12, height=3)
  
}

