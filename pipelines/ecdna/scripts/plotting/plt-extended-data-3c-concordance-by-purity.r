## Plot ecDNA detection concordance by binned tumor purity 
libs = c('optparse', 'ggplot2', 'ggpubr')
invisible(suppressPackageStartupMessages(sapply(libs, require, character.only=T)))
options(width=200, scipen=999)



## Get arguments
option_list = list(
  make_option(c("-a", "--amplicons"),       type='character', help="Merged amplicon summary"),
  make_option(c("-T", "--tm_file"),         type='character', help="Tumor-model pairs file"),
  make_option(c("-t", "--tn_file"),         type='character', help="Tumor-normal pairs file"),
  make_option(c("-p", "--purity_ploidy"),   type='character', help="Consensus purity/ploidy values"),
  make_option(c("-m", "--metadata"),        type='character', help="Patient metadata"),
  make_option(c("-o", "--out_file"),        type='character', help="Output SVG"))
opt = parse_args(OptionParser(option_list=option_list))


## Read tumor normal pairs 
tn = read.csv(opt$tn_file, h=T, stringsAsFactors=F)
tn$tumor = tn$tumor_aliquot_barcode
tn$normal = tn$normal_aliquot_barcode
tn = tn[, c('tumor', 'normal', 'pairId', 'participantId')]

## Read tumor-model pairs 
mt = read.csv(opt$tm_file, h=T, stringsAsFactors=F, sep='\t')
colnames(mt) = tolower(colnames(mt))

## Harmonize sample list so we're only using paired tumors/models
sample.list = intersect(unlist(mt[,c('tumor', 'model')]), tn$tumor)
tn = tn[tn$tumor %in% sample.list, ]
mt = mt[mt$tumor %in% sample.list | mt$model %in% sample.list, ]

## Read amplicons, filter for ecDNA in sample list 
aa = read.csv(opt$amplicons, h=T, stringsAsFactors=F, sep='\t')
aa = aa[aa$Classification == 'ecDNA', ]
aa$Sample.name = gsub('--.*', '', aa$Sample.name)
aa = aa[aa$Sample.name %in% tn$tumor, ]
aa$Case_ID = tn$participantId[match(aa$Sample.name, tn$tumor)]

## Read purity
pp = read.csv(opt$purity_ploidy, h=F, stringsAsFactors=F, col.names=c('tumor', 'purity','ploidy'))

## For each model-tumor comparison, annotate with tumor purity 
mt$tumor_purity = pp$purity[match(mt$tumor, pp$tumor)]
mt$model_purity = pp$purity[match(mt$model, pp$tumor)]

mt = data.frame(Sample.name=c(mt$tumor, mt$model), 
                tumor_purity=rep(mt$tumor_purity, 2),
                model_purity=rep(mt$model_purity, 2))

## Annotate AA with purity/ploidy
aa = merge(x=aa, y=mt, all.x=T, all.y=F, by='Sample.name')


## Bin purities
aa$tumor_purity_binned = cut(aa$tumor_purity, seq(0, 1, by=0.1))
aa$model_purity_binned = cut(aa$model_purity, seq(0, 1, by=0.1))

levels(aa$tumor_purity_binned) = gsub(',', '-', gsub('\\(|\\[|\\)|\\]', '', levels(aa$tumor_purity_binned)))
levels(aa$model_purity_binned) = gsub(',', '-', gsub('\\(|\\[|\\)|\\]', '', levels(aa$model_purity_binned)))

## Don't double-count concordant comparisons
aa = aa[is.na(aa$comparison_id) | !duplicated(aa$comparison_id), ]


## Format for plotting 
aa$concordance = ifelse(aa$model, 'Model only', 'Tumor only')
aa$concordance[aa$concordant_ecdna] = 'Intersecting'



##########
## Plot ##
##########

col = c(`Tumor only`='#ED4D44', 
        `Intersecting`='#8FBEDB',
        `Model only`='#D7A8B2')

aa$concordance = factor(aa$concordance, levels=rev(names(col)))


svg(opt$out_file, width=8, height=8)

# Plot by binned tumor purity
ggplot(aa, aes(x=tumor_purity_binned, fill=concordance)) +
              geom_bar(stat='count', position='stack') +
              scale_fill_manual(values=col) + 
              scale_x_discrete(drop=FALSE, name='Tumor purity bin window') +
              ylab('Number of ecDNAs detected') +
              theme_bw() + 
              guides(fill=guide_legend(nrow=3, reverse=T)) +
              theme(text=element_text(size=12, family='Arial Unicode MS'),
                    panel.grid.major=element_blank(),
                    panel.grid.minor=element_blank(),
                    strip.background=element_blank(),
                    panel.spacing=unit(5, 'mm'), 
                    plot.margin=unit(c(5.5, 8, 5.5, 8), "pt"),
                    legend.title=element_blank(),
                    legend.position='right')

dev.off()
message(opt$out_file)
