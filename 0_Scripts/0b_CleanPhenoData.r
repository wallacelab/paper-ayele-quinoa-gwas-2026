#! /usr/bin/Rscript

# Clean phenotype files and do QC analysis
library(argparse)
library(tidyverse)
library(performance)
library(emmeans)
library(ggpubr)

# Arguments
parser=ArgumentParser()
parser$add_argument("-i", "--infile", help="Input file of standardized pheno data (csv)")
parser$add_argument("-o", "--outprefix", help="Output file prefix")
parser$add_argument("-p", "--problem-traits", nargs="*", help="Traits that show a specific problem")
parser$add_argument("-m", "--prob-min-row", type="integer", default=0, help="Min row to filter out from Problem traits")
parser$add_argument("-x", "--prob-max-row", type="integer", default=0, help="Max row to filter out from Problem traits")
parser$add_argument("-q", "--iqr-cutoff", type="double", default=1.7, help="Cutoff for points to be considered an outlier by IQR")
parser$add_argument("-z", "--zscore-cutoff", type="double", default=3, help="Cutoff for points to be considered an outlier by Zscore")
args=parser$parse_args()
# setwd('/home/jgwall/Projects/Misc/AddissuAyele_QuinoaGBS/0_Phenotypes/')  # Debug working directory
# args=parser$parse_args(c("-i","0a_phenos.2024.csv",
#                          "-p", "PhiPS52", "ETR52", "VPD47",
#                          "-m", "1", "-x", "10",
#                          "-o",'99_tmp'))


# Load data
cat("Doing quality control on data from",args$infile,"\n")
raw = read.csv(args$infile) %>%
  mutate(Taxa = as.character(Taxa))

# Get list of traits
metacols = c("Taxa", "TRT") # Metadata columns
traits = names(raw) %>% 
  setdiff(metacols) %>%
  sort()

# Alphabetize traits (makes diagnostics easier to compare)
adjusted = raw %>%
  relocate(all_of(metacols), all_of(traits))

################
# Basic Cleaning
################


cat("Checking for outliers\n")

# Clean traits where first bit of TRT1 always shows a problem
if(length(args$problem_traits)>0){
  cat("Removing problem values for", args$problem_traits, "in rows",
      args$prob_min_row, "to", args$prob_max_row, "\n")
  problem_rows = args$prob_min_row:args$prob_max_row
  adjusted[problem_rows, args$problem_traits] = NA
}

# Use the performance package to check for outliers; since detects row-wise, split by trait
longdata = adjusted %>%
  pivot_longer(cols=all_of(traits), names_to="trait", values_to="value") %>%
  mutate(set = paste(trait, TRT, sep="|"))
# Clean via outlier detection
cat("\tFiltering outliers based on IQR cutoff",args$iqr_cutoff,"and z-score cutoff", args$zscore_cutoff, "\n")
sets = unique(longdata$set)
cleaned = lapply(sets, function(s){
  mydata = longdata %>% 
    filter(set==s) %>%
    filter(!is.na(value))  # Remove missing values from above; will be restored on pivot_wider()
  outliers = check_outliers(mydata, method=c("zscore", "iqr", "mcd"),
                            threshold = list(iqr=args$iqr_cutoff, zscore=args$zscore_cutoff)) %>%
    as.logical()
  mydata$value[outliers] = NA
  return(mydata)
})
filtered = cleaned %>%
  bind_rows() %>%
  select(-set) %>%
  pivot_wider(id_cols = all_of(metacols), names_from="trait", values_from="value")

#################
# Combine across reps
#################

# Make BLUEs
models = lapply(traits, function(mytrait){
  myformula = paste(mytrait, "~ TRT + Taxa")
  lm(myformula, data=filtered)
})
names(models) = traits

# Extract BLUEs with EMmeans (probably overkill)
taxa = unique(raw$Taxa) %>% sort()
blues = lapply(traits, function(mytrait){
  mymodel = models[[mytrait]]
  mymeans = emmeans(mymodel, specs="Taxa")
  
  # Return data frame
  myblues = mymeans %>%
    as.data.frame() %>%
    select(Taxa, emmean)
  names(myblues)[2] = mytrait
  return(myblues)
})
blues = blues %>%
  reduce(full_join, by="Taxa") # Combine into a single data frame

# Output BLUEs in csv format (wait on GAPIT format until after boxcox trans, etc.)
write.csv(filtered, file=paste(args$outprefix, ".filtered.csv", sep=""), row.names=FALSE)
write.csv(blues, file=paste(args$outprefix, ".blues.csv", sep=""), row.names=FALSE)
# blues.gapit = blues %>% rename("<Taxa>" = "Taxa")
# write.table(blues.gapit, file=paste(args$outprefix, ".gapit.txt", sep=""), row.names=FALSE,
#             col.names=TRUE, sep='\t', quote=FALSE)



##############
# Diagnostic plots
##############

cat("Outputting diagnostic plots (and flagging any problems in TRT1)\n")

# Make tidy data for plotting
raw.plotdata = raw %>%
  pivot_longer(cols=all_of(traits), names_to="Trait", values_to="Value")
filt.plotdata = filtered %>%
  pivot_longer(cols=all_of(traits), names_to="Trait", values_to="Value")
blues.plotdata = blues %>%
  pivot_longer(cols=all_of(traits), names_to="Trait", values_to="Value")

# Function to make all diagnostic plots for a specific trait
make_plots = function(mytrait){
  myraw = raw.plotdata %>% 
    filter(Trait==mytrait)
  myfilt = filt.plotdata %>% 
    filter(Trait==mytrait)
  myblues = blues.plotdata %>%
    filter(Trait==mytrait)
  
  # Find outliers
  outliers = subset(myraw, is.na(myfilt$Value))
  # Find problematic ones in TRT 1
  if(mytrait %in% args$problem_traits){
    problems = subset(outliers, outliers$Taxa %in% raw$Taxa[problem_rows] &
                        outliers$TRT=="TRT1")
  }else{
    problems=outliers[FALSE,] # 0-row data frame
  }
  
  # Raw data plot
  scatter.raw = ggplot(myraw) +
    aes(x=Taxa, y=Value, color=TRT) +
    geom_point(color="red", shape=1, data=outliers, size=4, stroke=1.8) + 
    geom_point(color="black", shape=1, data=problems, size=4, stroke=2) + # Paint overtop
    geom_point(alpha=0.5) +
    theme(axis.text.x = element_blank(), title = element_text(face="bold"),
          legend.position="none") +
    labs(title=paste(mytrait, "(raw)"))
  
  # Raw histogram
  hist.raw = ggplot(myraw) +
    aes(x=Value, fill=TRT) +
    geom_histogram(bins=30) +
    facet_wrap(~ TRT, nrow=1) +
    labs(title="Raw data distribution by block")
  
  # Filtered data plot
  scatter.filt = ggplot(myfilt) +
    aes(x=Taxa, y=Value, color=TRT) +
    geom_point(alpha=0.5) +
    theme(axis.text.x = element_blank(), title = element_text(face="bold"),
          legend.position="none") +
    labs(title=paste(mytrait, "(filtered)"))
  
  # Filtered histogram
  hist.filt = ggplot(myfilt) +
    aes(x=Value, fill=TRT) +
    geom_histogram(bins=30) +
    facet_wrap(~ TRT, nrow=1) +
    labs(title="Filtered data distribution by block")
  
  # BLUEs scatterplot
  scatter.blues = ggplot(myblues) +
    aes(x=Taxa, y=Value) +
    geom_point(alpha=0.5, color="black") +
    theme(axis.text.x = element_blank()) +
    labs(title="BLUEs (from filtered)")
  
  # BLUEs histogram
  hist.blues = ggplot(myblues) +
    aes(x=Value) +
    geom_histogram(color="black", bins=30) +
    labs(title="BLUEs distribution (from filtered)")
  
  # # GVLMA tests - Saving for data cleaning step
  # mymodel = models[[mytrait]]
  # mytest = gvlma(mymodel) %>% 
  #   summary() %>%
  #   rownames_to_column("test") %>%
  #   mutate(Decision = sub(Decision, pattern="Assumptions acceptable.", repl="Yes")) %>%
  #   mutate(Decision = sub(Decision, pattern="Assumptions NOT satisfied!", repl="NO")) %>%
  #   mutate(p=round(`p-value`, digits=2)) %>%
  #   mutate(yval=1:n())
  # testplot = ggplot(mytest) + 
  #   aes(y=yval, color=Decision) +
  #   theme_void() +
  #   theme(legend.position="none") + 
  #   geom_text(mapping=aes(label=test), x=1, hjust=1) +
  #   geom_text(mapping=aes(label=p), x=1.25, hjust=0) +
  #   geom_text(mapping=aes(label=Decision), x = 1.75) +
  #   lims(x=c(-1, 2.5), y=c(0,6)) +
  #   scale_color_manual(values=c("Yes"="darkblue", "NO"= "darkred")) +
  #   labs(title="GVLMA tests (with p-value and decision)")
  
  # Arrange together
  #myplots = ggarrange(scatter.raw, hist.raw, scatter.blues, hist.blues, testplot,
  myplots = ggarrange(scatter.raw, hist.raw, scatter.filt, hist.filt, scatter.blues, hist.blues,
                      nrow=1, widths=c(2,3,2,3,2,2))
  return(myplots)
}




# Combine and output
allplots = lapply(traits, make_plots)
#allplots = lapply(c("DBI", "ETR40", "HSW"), make_plots) # For debugging
combined = ggarrange(plotlist=allplots, ncol=1) # 1 column so all get stacked on top of each other
ggsave(combined, file=paste(args$outprefix, ".diagnistics.png", sep=""),
       width=20, height=3 * length(allplots), limitsize = FALSE, dpi=200)

