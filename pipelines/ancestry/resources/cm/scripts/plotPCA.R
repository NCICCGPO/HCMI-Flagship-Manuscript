#!/usr/bin/Rscript
# r script to plot pca results for a set of references + one study sample

args <- commandArgs(TRUE)
inputFile <- args[1]
admixRes <- args[2]
outputFile <- args[3]

library(ggplot2)
library(RColorBrewer)
library(gridExtra)
library(grid)


grid_arrange_shared_legend <- function(...) {
    plots <- list(...)
    g <- ggplotGrob(plots[[1]] + theme(legend.position="bottom"))$grobs
    legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
    lheight <- sum(legend$height)
    grid.arrange(
        do.call(arrangeGrob, lapply(plots, function(x)
            x + theme(legend.position="none"))),
        legend,
        ncol = 1,
        heights = unit.c(unit(1, "npc") - lheight, lheight))
}

myColors <- brewer.pal(5,"Set1")
myColors<-c(myColors,"gray")
names(myColors) <- c("EAS","SAS","EUR","AMR","AFR","study sample")
colScale <- scale_colour_manual(name = "population",values = myColors)

pca <- read.table(inputFile,header=F,as.is=T)
res <- read.table(admixRes,header=F,as.is=T)

res <- res[,c(1,3)]
colnames(res) <- c("pop","sampleID")
pca <- merge(pca,res,by.x="V2",by.y="sampleID",all.x=TRUE)

pca$study <- FALSE
pca$study[is.element(pca$pop,"-")] <- TRUE
pca$pop[is.element(pca$pop,"-")] <- "study sample"



p1<-ggplot(pca,aes(x=V3,y=V4,colour=pop)) + geom_point() + theme_bw() + xlab("EV1") + ylab("EV2") + colScale +
geom_point(data=pca[pca$study,],aes(x=V3,y=V4),color="black",shape=3,size=3)

p2<-ggplot(pca,aes(x=V4,y=V5,colour=pop)) + geom_point() + theme_bw() + xlab("EV2") + ylab("EV3") + colScale +
geom_point(data=pca[pca$study,],aes(x=V4,y=V5),color="black",shape=3,size=3)

p3<-ggplot(pca,aes(x=V5,y=V6,colour=pop)) + geom_point() + theme_bw() + xlab("EV3") + ylab("EV4") + colScale +
geom_point(data=pca[pca$study,],aes(x=V5,y=V6),color="black",shape=3,size=3)

p4<-ggplot(pca,aes(x=V6,y=V7,colour=pop)) + geom_point() + theme_bw() + xlab("EV4") + ylab("EV5") + colScale +
geom_point(data=pca[pca$study,],aes(x=V5,y=V6),color="black",shape=3,size=3)

pdf(paste0(outputFile), width=8, height=8)
grid_arrange_shared_legend(p1, p2, p3, p4)
dev.off()

