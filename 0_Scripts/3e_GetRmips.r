#! /usr/bin/Rscript

# Compile permutations and resamples to get RMIPs
library(argparse)
library(ggplot2)
library(dplyr)
library(parallel)
options(stringsAsFactors=F)

# Arguments
parser=ArgumentParser()
parser$add_argument("-p", "--perms", nargs="*", help="List of input files of GAPIT results (permutations)")
parser$add_argument("-r", "--resamples", nargs="*", help="List of input files of GAPIT results (resample iterations)")
parser$add_argument("--perm-p-cutoff", type='double', default=0.01, help="P-value cutoff (determined from permutations) to be considered a hit in a resample run")
parser$add_argument("--rmip-cutoff", type='double', default=0.05, help="RMIP cutoff to be consider a genuine hit")
parser$add_argument("-o", "--outprefix", help="Output file prefix")
parser$add_argument("-t", "--trait", default="unnamed trait", help="Trait name")
parser$add_argument("-n", "--num-cores", default=1, help="Number of parallel cores to use to load data")
args=parser$parse_args()
# setwd('/home/jgwall/Projects/Misc/AddissuAyele_QuinoaGBS/3_GWAS/')
# args=parser$parse_args(c("-p","3b_Permutations.FarmCPU/PhiPS56/perm1/GAPIT.Association.GWAS_Results.FarmCPU.PhiPS56.csv.gz",
#                        "3b_Permutations.FarmCPU/PhiPS56/perm2/GAPIT.Association.GWAS_Results.FarmCPU.PhiPS56.csv.gz",
#                        "-r", "3d_Resampling.FarmCPU/PhiPS56/resample1/GAPIT.Association.GWAS_Results.FarmCPU.PhiPS56.csv.gz",
#                        "3d_Resampling.FarmCPU/PhiPS56/resample2/GAPIT.Association.GWAS_Results.FarmCPU.PhiPS56.csv.gz",
#                        "-o", "99_tmp", "-n", "2"))


cat("Compiling permutations and resamplings for",args$outprefix,"on",args$num_cores,"cores\n")

# Load permutations
cat("\tLoading permutations from",length(args$perms),"input files\n")
perms = mclapply(args$perms, mc.cores=args$num_cores, FUN=function(infile){
  mygwas = read.csv(infile, colClasses="character") # Read as character to speed up
  mygwas$P.value = as.numeric(mygwas$P.value)
  best = which.min(mygwas$P.value)
  return(mygwas[best,])
})
perms = sort(do.call(rbind, perms)$P.value)

# Load resamples; process at this stage so as to not save a bunch of unnecessary data in memory
cat("\tLoading resamples from",length(args$resamples),"input files\n")
resamples = mclapply(args$resamples, mc.cores=args$num_cores, FUN=function(infile){
    mydata = read.csv(infile, colClasses="character")
    mydata$P.value = as.numeric(mydata$P.value)
    
    # # Get empirical p-value
    # mydata$empirical_p = sapply(mydata$P.value, function(p){
    #     sum(p>=perms)/length(perms)
    # })
    
    # Below is equivalent to above; need to determine which is faster
    mydata$empirical_p = findInterval(mydata$P.value, vec=perms) / length(perms)
    
    
    # Filter based on input arguments
    mydata=subset(mydata, mydata$empirical_p <= args$perm_p_cutoff)
    return(mydata[,c('SNP', 'Chr', 'Pos', 'empirical_p')])
    #return(mydata)
})

# Get RMIPs
cat("\tCalculating RMIPs\n")
resamples.all = do.call(rbind, resamples)
rmip = paste(resamples.all$SNP, resamples.all$Chr, resamples.all$Pos, sep="###")
rmip = table(rmip) / length(args$resamples)
cat("\t\tTotal of", sum(rmip >= args$rmip_cutoff), "traits pass RMIP cutoff\n")

# Separate for plotting (and make null results if none present)
if(length(rmip)==0){
    toplot = data.frame(snp=character(), chr=numeric(), pos=numeric(), rmip=numeric())
    myplot= ggplot() + 
      annotate("text", x = 4, y = 25, size=8,
               label = paste("No hits for",args$trait), color="red") + 
      theme_void()
}else{
  # Make data table
  toplot = strsplit(names(rmip), split="###")
  toplot = data.frame(do.call(rbind, toplot))
  names(toplot) = c("snp", "chr", "pos")
  toplot$chr = as.numeric(toplot$chr)
  toplot$pos = as.numeric(toplot$pos)
  toplot$rmip = as.numeric(rmip)
  
  # Graph results
  toplot$significant = ifelse(toplot$rmip >= args$rmip_cutoff, yes="yes", no="no")
  myplot = ggplot(toplot, mapping=aes(x=pos, y=rmip, color=significant)) +
    geom_point() + 
    facet_wrap(~chr, nrow=2) +
    ggtitle(args$trait) +
    scale_color_manual(values=c("yes"="red", "no"="black"))
}
    

# Output 
cat("\tWriting out results\n")
toplot = toplot[order(toplot$rmip, toplot$chr, toplot$pos),]
write.csv(toplot, file=paste(args$outprefix, ".csv", sep=""), row.names=F)
ggsave(myplot, file=paste(args$outprefix, ".png", sep=""), width=10, height=6, dpi=300)
