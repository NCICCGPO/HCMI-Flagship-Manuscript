# script to parse .traw	file of	genotypes to PLINK format

args <-	  commandArgs(TRUE)
(fname <- args[1])
(outFname <- args[2])

library(readr)
library(GWASTools)
library(SNPRelate)

genos <- read_delim(fname,col_names=TRUE,delim="\t")

scanAnnot <- data.frame(scanID=colnames(genos)[-c(1:5)],
	  subjectID=colnames(genos)[-c(1:5)])
	  
colnames(scanAnnot) <- c("scanID","subjectID")

scanAnnot$intID <- 1:nrow(scanAnnot)
scanAnnot$scanID <- scanAnnot$intID

numGenos <- as.matrix(genos[,-c(1:5)])
numGenos <- as.integer(numGenos)

gMat <- MatrixGenotypeReader(genotype=matrix(numGenos,nrow=nrow(genos),ncol=ncol(genos)-5,byrow=FALSE),snpID=1:nrow(genos),
chromosome=as.integer(genos$chr),position=as.integer(genos$pos),scanID=scanAnnot$scanID)

snpAnnot <- data.frame(snpID=1:nrow(genos),rsID=as.character(genos$rsID),chromosome=as.integer(genos$chr),position=as.integer(genos$pos),
alleleA=as.character(genos$alleleA),alleleB=as.character(genos$alleleB))
snpAnnot$alleleA <- as.character(snpAnnot$alleleA)
snpAnnot$alleleB <- as.character(snpAnnot$alleleB)

gData <- GenotypeData(gMat,scanAnnot=ScanAnnotationDataFrame(scanAnnot),snpAnnot=SnpAnnotationDataFrame(snpAnnot))
vcfWrite(gData,id.col="rsID",vcf.file=outFname,sample.col="subjectID")
