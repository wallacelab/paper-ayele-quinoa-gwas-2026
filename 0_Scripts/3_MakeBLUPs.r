#! /usr/bin/Rscript

# Check and unify genotype and phenotype files
library(argparse)
library(tidyverse)
library(emmeans)
library(gvlma)
library(ggpubr)

# Arguments
parser=ArgumentParser()
parser$add_argument("-i", "--infile", help="Input file of phenotype data (csv)")
parser$add_argument("-o", "--outprefix", help="Output file prefix")
args=parser$parse_args()
# setwd('/home/jgwall/Projects/Misc/AddissuAyele_QuinoaGBS/')  # Debug working directory
# args=parser$parse_args(c("-i","0_Scripts/2024_Greenhouse data.csv",
#                          "-o",'99_tmp'))


# Load data
cat("Loading raw data to convert to BLUEs\n")
raw=read.csv(args$infile) %>%
  mutate(Taxa = as.character(Taxa), ENTRY=as.character(ENTRY), TRT=as.character(TRT)) # Make characters so treated as factors
traits = names(raw) %>% setdiff(c("Taxa", "ENTRY", "TRT")) # Get list of traits

# Make BLUEs
models = lapply(traits, function(mytrait){
  myformula = paste(mytrait, "~ TRT + Taxa")
  lm(myformula, data=raw)
})
names(models) = traits

# Extract BLUEs with EMmeans
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

# Output BLUEs in csv and GAPIT formats
write.csv(blues, file=paste(args$outprefix, ".csv", sep=""), row.names=FALSE)
blues.gapit = blues %>% rename("<Taxa>" = "Taxa")
write.table(blues.gapit, file=paste(args$outprefix, ".gapit.txt", sep=""), row.names=FALSE,
            col.names=TRUE, sep='\t', quote=FALSE)


##############
# Diagnostic plots
##############

# Make tidy data for plotting
raw.plotdata = raw %>%
  pivot_longer(cols=all_of(traits), names_to="Trait", values_to="Value")
blues.plotdata = blues %>%
  pivot_longer(cols=all_of(traits), names_to="Trait", values_to="Value")

# Function to make all diagnostic plots for a specific trait
make_plots = function(mytrait){
  myraw = raw.plotdata %>% 
    filter(Trait==mytrait)
  myblues = blues.plotdata %>%
    filter(Trait==mytrait)
  
  # Raw data plot
  scatter.raw = ggplot(myraw) +
    aes(x=Taxa, y=Value, color=TRT) +
    geom_point(alpha=0.5) +
    theme(axis.text.x = element_blank(), title = element_text(face="bold")) +
    labs(title=mytrait)
  
  # Raw histogram
  hist.raw = ggplot(myraw) +
    aes(x=Value, fill=TRT) +
    geom_histogram() +
    facet_wrap(~ TRT, nrow=1) +
    labs(title="Raw data distribution by block")
  
  # BLUEs scatterplot
  scatter.blues = ggplot(myblues) +
    aes(x=Taxa, y=Value) +
    geom_point(alpha=0.5, color="black") +
    theme(axis.text.x = element_blank()) +
    labs(title="BLUEs")
  
  # BLUEs histogram
  hist.blues = ggplot(myblues) +
    aes(x=Value) +
    geom_histogram(color="black") +
    labs(title="BLUEs distribution")
  
  # GVLMA tests
  mymodel = models[[mytrait]]
  mytest = gvlma(mymodel) %>% 
    summary() %>%
    rownames_to_column("test") %>%
    mutate(Decision = sub(Decision, pattern="Assumptions acceptable.", repl="Yes")) %>%
    mutate(Decision = sub(Decision, pattern="Assumptions NOT satisfied!", repl="NO")) %>%
    mutate(p=round(`p-value`, digits=2)) %>%
    mutate(yval=1:n())
  testplot = ggplot(mytest) + 
    aes(y=yval, color=Decision) +
    theme_void() +
    theme(legend.position="none") + 
    geom_text(mapping=aes(label=test), x=1, hjust=1) +
    geom_text(mapping=aes(label=p), x=1.25, hjust=0) +
    geom_text(mapping=aes(label=Decision), x = 1.75) +
    lims(x=c(-1, 2.5), y=c(0,6)) +
    scale_color_manual(values=c("Yes"="darkblue", "NO"= "darkred")) +
    labs(title="GVLMA tests (with p-value and decision)")
    
  # Arrange together
  myplots = ggarrange(scatter.raw, hist.raw, scatter.blues, hist.blues, testplot,
                      nrow=1, widths=c(2,3,2,2,2))
  return(myplots)
}




# Combine and output
allplots = lapply(traits, make_plots)
#allplots = lapply(c("DBI", "ETR40", "HSW"), make_plots) # For debugging
combined = ggarrange(plotlist=allplots, ncol=1) # 1 column so all get stacked on top of each other
ggsave(combined, file=paste(args$outprefix, ".diagnistics.png", sep=""),
       width=20, height=3 * length(allplots), limitsize = FALSE, dpi=200)


