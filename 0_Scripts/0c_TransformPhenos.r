#! /usr/bin/Rscript

# Transform supplied phenos with Boxcox transformation to find best fit

# Load libraries
library(argparse)
library(tidyverse)
library(MASS)
library(ggpubr)
library(scales)

# Arguments
parser=ArgumentParser()
parser$add_argument("-i", "--infile", help="Input file of standardized pheno data (csv)")
parser$add_argument("-o", "--outprefix", help="Output file prefix")
args=parser$parse_args()
# setwd('/home/jgwall/Projects/Misc/AddissuAyele_QuinoaGBS/0_Phenotypes/')  # Debug working directory
# args=parser$parse_args(c("-i","0b_phenos.2025.filtered.csv",
#                          "-o",'99_tmp'))

cat("Finding best transformation of data in", args$infile,"by box-cox\n")
cat("\tNote: All values will be rescaled to range of 1-100 to make more effective\n")
phenos = read.csv(args$infile)

# Get list of traits
metacols = c("Taxa", "TRT") # Metadata columns
traits = names(phenos) %>% 
  setdiff(metacols) %>%
  sort()

# Make data long and split into subsets based on trait and treatment
longdata = phenos %>%
  pivot_longer(cols=all_of(traits), names_to="trait", values_to="value")
subsets =split(longdata, f=list(longdata$TRT, longdata$trait))

# Fit linear model before and after boxcox transformation to each variable
transformed = list()
for(i in 1:length(subsets)){ # Have to do a for loop b/c boxcox issue
  # Make sure response data is positive
  mydata = subsets[[i]]
  mydata = subset(mydata, !is.na(mydata$value)) # Remove missing data
  mydata$value_trans = rescale(mydata$value, to=c(1,100)) # Rescale so values are between 1 and 100 so boxcox works well
  mymodel = lm("value_trans ~ 1", data = mydata)
  
  # Determine best boxcox transformation
  bc = boxcox(mymodel, lambda=seq(-4, 4, 0.05), plotit = FALSE)
  best = which.max(bc$y)
  lambda = bc$x[best]
  
  # Transform data
  if(lambda == 0){
    mydata$value_trans = log(mydata$value_trans)
  }else{
    mydata$value_trans = (mydata$value_trans ^ lambda - 1)/lambda
  }
  mydata$lambda = lambda
  transformed[[i]] = mydata
}
transformed = bind_rows(transformed)

#########
# Diagnostic histograms before and afterward
#########

cat("Making diagnostic plots\n")

# Function to make all diagnostic plots for a specific trait
make_plots = function(mytrait){
  
  mydata = subset(transformed, trait==mytrait)
  
  # Raw histogram
  hist.orig = ggplot(mydata) +
    aes(x=value, fill=TRT) +
    geom_histogram(bins=30) +
    facet_wrap(~ TRT, nrow=1, scales="free_x") +
    labs(title=paste(mytrait, " - Original data"))
  
  # Filtered histogram
  lambdas = mydata %>%
    dplyr::select(TRT, lambda) %>%
    unique() %>%
    mutate(label = paste("λ =", lambda))
  hist.trans = ggplot(mydata) +
    aes(x=value_trans, fill=TRT) +
    geom_histogram(bins=30) +
    geom_text(data=lambdas, mapping=aes(label=label), x=-Inf, y=Inf,
              hjust=-0.1, vjust=1.5) +
    facet_wrap(~ TRT, nrow=1, scales="free_x") +
    labs(title=paste(mytrait, " - Boxcox-transformed"))
  
  # Arrange together
  myplots = ggarrange(hist.orig, hist.trans)
  return(myplots)
}

# Combine and output
allplots = lapply(traits, make_plots)
combined = ggarrange(plotlist=allplots, ncol=1) # 1 column so all get stacked on top of each other
ggsave(combined, file=paste(args$outprefix, ".trans_plots.png", sep=""),
       width=8, height=2 * length(allplots), limitsize = FALSE, dpi=200)


##############
# Output transformed data
###############

outdata = transformed %>%
  dplyr::select(-value, -lambda) %>%
  pivot_wider(id_cols = all_of(metacols), names_from="trait", values_from="value_trans")
write.csv(outdata, file=paste(args$outprefix, ".transformed.csv", sep=""), row.names=FALSE)

lamdbavals = transformed %>%
  dplyr::select(trait, TRT, lambda) %>%
  unique() %>%
  mutate(lambda = round(lambda, digits=3)) # Taking care of a floating point display error
write.csv(lamdbavals, file=paste(args$outprefix, ".transformed.lambdas.csv", sep=""), row.names=FALSE)

