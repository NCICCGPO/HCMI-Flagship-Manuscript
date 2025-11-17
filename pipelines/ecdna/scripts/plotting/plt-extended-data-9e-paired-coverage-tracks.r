## Plot coverage for the union of two concordant tumor/model ecDNAs 
libs = c('optparse', 'GenomicRanges', 'gGnome', 'gTrack', 'gUtils', 'rtracklayer')
invisible(suppressPackageStartupMessages(sapply(libs, require, character.only=T)))



## Read merged tumor-model feature file and take the union of the intervals 
read.feature.file = function(f, amps) {

  x = read.csv(f, h=F, stringsAsFactors=F, sep='\t', col.names=c('intervals', 'graph'))

  x = x[gsub('_intervals.bed', '', basename(x$intervals)) %in% amps, ]

  res = Reduce(c, lapply(x$intervals, rtracklayer::import))

  return(res)

}



## Alpha function borrowed from mski lab at some point
alpha = function(col, alpha) {
  col.rgb = col2rgb(col)
  out = rgb(red = col.rgb['red', ]/255, green = col.rgb['green', ]/255, blue = col.rgb['blue', ]/255, alpha = alpha)
  names(out) = names(col)
  return(out)
}



## Get arguments
option_list = list(
  make_option(c("-T", "--tumor"),          type='character', help="Tumor name"),
  make_option(c("-M", "--model"),          type='character', help="Tumor name"),
  make_option(c("-a", "--amp1"),           type='character', help="ID of amplicon 1"),
  make_option(c("-A", "--amp2"),           type='character', help="ID of amplicon 2"),
  make_option(c("-c", "--tumor_coverage"), type='character', help="Tumor coverage RDS"),
  make_option(c("-C", "--model_coverage"), type='character', help="Model coverage RDS"),
  make_option(c("-F", "--feat2graph"),     type='character', help="Features to graph table for this tumor-model comparison"),
  make_option(c("-f", "--field"),          type='character', help="Field in --coverage file to use when plotting coverage", default='foreground'),
  make_option(c("-g", "--genes"),          type='character', help="Gene list for gencode track"),
  make_option(c("-o", "--out_file"),       type='character', help="Figure output file (SVG)"))
opt = parse_args(OptionParser(option_list=option_list))



## Init gene track 
# GENCODE = '/gpfs/commons/projects/nepc/analysis/jabba-latest/annotations/all-genes/gencode.composite.collapsed.rds'


if (!is.null(opt$genes) & file.exists(opt$genes)) {
  gencode = track.gencode(gencode=GENCODE, 
                          cached.dir=dirname(GENCODE), 
                          cached.path=GENCODE, 
                          cex.label=1, 
                          xaxis.cex.label=1, 
                          xaxis.unit=1e6, 
                          xaxis.suffix='MB',
                          genes=scan(opt$genes, what=character())
                          )
} else if (!is.null(opt$genes) & !file.exists(opt$genes)) {
  gencode = track.gencode(gencode=GENCODE, 
                          cached.dir=dirname(GENCODE), 
                          cached.path=GENCODE, 
                          cex.label=1.2, 
                          xaxis.cex.label=1, 
                          xaxis.unit=1e6, 
                          xaxis.suffix='MB',
                          genes=unlist(strsplit(opt$genes, ',')))
} else {
  gencode = track.gencode(gencode=GENCODE, 
                          cached.dir=dirname(GENCODE), 
                          cached.path=GENCODE, 
                          cex.label=1.2, 
                          xaxis.cex.label=1, 
                          xaxis.unit=1e6, 
                          xaxis.suffix='MB')
}
gencode$height = 4
gencode$legend = FALSE



## Read feature file and take union of ecDNA intervals 
bed = read.feature.file(f=opt$feat2graph, amps=c(opt$amp1, opt$amp2))



## Read coverage
tumor.cov = readRDS(opt$tumor_coverage)
seqlevelsStyle(tumor.cov) = 'UCSC'
tumor.cov.gt = gTrack(tumor.cov, 
                      y.field=opt$field, 
                      col=alpha('black', 0.2), 
                      lwd.border=1, 
                      name='Tumor', 
                      ylab=' ')

model.cov = readRDS(opt$model_coverage)
seqlevelsStyle(model.cov) = 'UCSC'
model.cov.gt = gTrack(model.cov, 
                      y.field=opt$field, 
                      col=alpha('black', 0.2), 
                      lwd.border=1, 
                      name='Model', 
                      ylab=' ')



## Build final gTrack
plt.gt = c(gencode, model.cov.gt, tumor.cov.gt)



##########
## Plot ## 
##########

svg(opt$out_file, width=2.53, height=2.53)

par(family='Arial Unicode MS', cex=0.45)
plot(plt.gt, 
     bed+1E5, 
     chr.sub=F,
     xaxis.interval=1E7,
     y.grid.lwd=0.4,
     sep.lwd=0.6,
     sep.lty=3)

dev.off()

message(opt$out_file)


