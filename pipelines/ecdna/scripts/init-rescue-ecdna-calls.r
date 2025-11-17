## Rescue filtered ecDNAs if they were part of a concordant tumor/model overlap
libs = c('optparse', 'reshape2')
invisible(suppressPackageStartupMessages(sapply(libs, require, character.only=T)))
options(width=200, scipen=999)



parse.genes = function(x) {

  x = gsub('\\]', ')', gsub('\\[', 'c(', x))
  x = lapply(x, function(y) eval(parse(text=y)))
  
  return(x)

}



## Get arguments
option_list = list(
  make_option(c("-i", "--in_file_amplicons"),        type='character', help="Unfiltered amplicon table"),
  make_option(c("-I", "--in_file_concordant_ecdna"), type='character', help="Concordant tumor/model ecdnas"),
  make_option(c("-o", "--out_file"),                 type='character', help="Output TSV"))
opt = parse_args(OptionParser(option_list=option_list))



## Read amplicons 
aa = read.csv(opt$in_file_amplicons, h=T, stringsAsFactors=F, sep='\t')

## Reformat genes for readability
aa$All.genes = sapply(parse.genes(aa$All.genes), paste, collapse=',')
aa$Oncogenes = sapply(parse.genes(aa$Oncogenes), paste, collapse=',')

## Reformat intervals for easier downstream parsing 
aa$Location = sapply(parse.genes(aa$Location), paste, collapse=',')

## Read concordant ecDNAs 
conc = read.csv(opt$in_file_concordant_ecdna, h=T, stringsAsFactors=F, sep='\t')

## Map between comparison IDs and amplicon IDs
id.map = melt(conc[, c('comparison_id', 'Amp1', 'Amp2')], id.var='comparison_id', variable.name='amp_12', value.name='amp_id')
aa$comparison_id = id.map$comparison_id[match(aa$Feature.ID, id.map$amp_id)]
aa$concordant_ecdna = !is.na(aa$comparison_id)


## Retain anything with no filters, OR an ecDNA that has a concordant call in its 
## paired tumor/model
aa = aa[aa$Filter.flag == 'None' | aa$concordant_ecdna, ]

## Remove columns with file paths 
aa = aa[, grep('file|JSON', colnames(aa), invert=T)]

## These fields aren't used 
aa = aa[, !colnames(aa) %in% c('Tissue.of.origin', 'Sample.type')]


## Write result 
write.table(aa, opt$out_file, row.names=F, col.names=T, quote=F, sep='\t')
message(opt$out_file)


## Detection summary post-rescue 
message('\n --- Post-rescue detection summary ---')
aa$model = factor(ifelse(aa$model, 'Model', 'Tumor'), levels=c('Tumor', 'Model'))

message('\nCall counts: ')
call.counts = as.data.frame.matrix(table(aa$Classification, aa$model))
write.table(call.counts, row.names=T, col.names=T, quote=F, sep='\t')

message('\nSample counts: ')
sample.counts = tapply(aa$Sample.name, list(aa$Classification, aa$model), function(x) length(unique(x)))
write.table(sample.counts, row.names=T, col.names=T, quote=F, sep='\t')

message('\necDNA count summary')
tapply(aa$Sample.name[aa$Classification=='ecDNA'], aa$model[aa$Classification=='ecDNA'], function(x) summary(as.numeric(table(x))))

message('\nConcordance summary')
table(aa$concordant_ecdna[aa$Classification=='ecDNA'], aa$model[aa$Classification=='ecDNA'])




message('\n-- Oncogene summary : TUMORS --')
oncogenes = strsplit(aa$Oncogenes[aa$Classification == 'ecDNA' & aa$model == 'Tumor'], ',')
message('Unique tumor oncogenes:', length(unique(unlist(oncogenes))))
write.table(head(as.data.frame(sort(table(unlist(oncogenes)), decreasing=T)), 25), sep='\t', row.names=F, quote=F)


message('\n-- Oncogene summary : MODELS --')
oncogenes = strsplit(aa$Oncogenes[aa$Classification == 'ecDNA' & aa$model == 'Model'], ',')
message('Unique model oncogenes:', length(unique(unlist(oncogenes))))
write.table(head(as.data.frame(sort(table(unlist(oncogenes)), decreasing=T)), 25), sep='\t', row.names=F, quote=F)


message('\n-- Oncogene summary : CONCORDANT PAIRS --')
oncogenes = strsplit(aa$Oncogenes[aa$Classification == 'ecDNA' & aa$concordant_ecdna], ',')
message('Unique model oncogenes:', length(unique(unlist(oncogenes))))
write.table(head(as.data.frame(sort(table(unlist(oncogenes))/2, decreasing=T)), 25), sep='\t', row.names=F, quote=F)


message('\n-- Oncogene summary : ALL --')
oncogenes = strsplit(aa$Oncogenes[aa$Classification == 'ecDNA' & aa$model == 'Model'], ',')
message('Unique oncogenes across all samples:', length(unique(unlist(oncogenes))))
