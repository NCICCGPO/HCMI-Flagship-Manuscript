
# R script to parse output from ADMIXTURE runs
# takes as argument one file name, the output from batch_runADMIXTURE.sh
# assumes:
# study ancestry is coded as -
# outputs just the study samples and their proportions
# file output is inputFileName.studyOnly

args <- commandArgs(TRUE)
resFile <- args[1]

res <- read.table(resFile,header=F,as.is=T)

for(i in 8:ncol(res)){
ancest <- res$V1[res[,i]>0.99][2]
colnames(res)[i] <- ancest
}

resSamp <- res[is.element(res$V1,"-"),c(3,8:ncol(res))]
colnames(resSamp)[1] <- "sampleID"

write.table(resSamp,file=paste0(resFile,".studyOnly"),quote=F,row.names=F)

