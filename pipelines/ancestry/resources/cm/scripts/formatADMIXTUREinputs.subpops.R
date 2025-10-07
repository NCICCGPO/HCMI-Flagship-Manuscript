#!/usr/bin/env Rscript
# script to format .pop input file for supervised ADMIXTURE analysis

args <- commandArgs(TRUE)
input <- args[1] # file path/root of .fam file for ADMIXTURE input
kgpops <- args[2]
popType <- args[3] # this is either "superpops" or one of AFR,EUR,SAS,EAS,AMR

fam <- read.table(paste0(input,".fam"), header=FALSE, as.is=TRUE)
pop <- fam$V1

popFile <- read.table(kgpops, header=TRUE, as.is=TRUE, sep="\t")

if(popType=="superpops"){
for(i in 1:nrow(popFile)){
  idx <- which(is.element(pop,popFile$Population.code[i]))
    pop[pop==popFile$Population.code[i]] <- popFile$Population.code[i]
    }
    pop[!is.element(pop,unique(popFile$Population.code))] <- "-"
    table(pop)
    }


write.table(pop,file=paste0(input,".pop"), row.names=FALSE, col.names=FALSE, quote=FALSE)

q("no")
    
