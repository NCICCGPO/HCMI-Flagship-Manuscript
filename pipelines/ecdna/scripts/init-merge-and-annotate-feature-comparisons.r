## Merge and annotate results of AmpliconClassifier's feature_similarity.py
libs = c('optparse')
invisible(suppressPackageStartupMessages(sapply(libs, require, character.only=T)))
options(width=200, scipen=999)



parse.genes = function(x) {

  x = gsub('\\]', ')', gsub('\\[', 'c(', x))
  x = lapply(x, function(y) eval(parse(text=y)))
  
  return(x)

}



## Extract complexity score for a given amplicon
get.complexity.score = function(f, amplicon_id) {

  x = read.csv(f, h=T, stringsAsFactors=F, sep='\t', check.names=F)
  res = x$`Complexity score`[x$`Feature ID` == amplicon_id]
  return(res)

}



## Get arguments
option_list = list(
  make_option(c("-T", "--tm_file"),                          type='character', help="Tumor-normal pairs file"),
  make_option(c("-A", "--aa_dir"),                           type='character', help="Top-level AmpliconArchitect output directory"),
  make_option(c("-t", "--tn_file"),                          type='character', help="Tumor-normal pairs file"),
  make_option(c("-a", "--amplicons"),                        type='character', help="Merged amplicon summary"),
  make_option(c("-m", "--min_jaccard_similarity"),           type='numeric',   help="Minimum Jaccard interval similarity to consider as concordant"),
  make_option(c("-f", "--fc_dir"),                           type='character', help="Dir with feature_similarity.py output"),
  make_option(c("-o", "--out_file"),                         type='character', help="Output TSV"))
opt = parse_args(OptionParser(option_list=option_list))



## Read tumor normal pairs 
tn = read.csv(opt$tn_file, h=T, stringsAsFactors=F)
tn$tumor = tn$tumor_aliquot_barcode
tn$normal = tn$normal_aliquot_barcode
tn$pair_name = tn$pairId


## Read tumor-model pairs 
mt = read.csv(opt$tm_file, h=T, stringsAsFactors=F, sep='\t')
colnames(mt) = tolower(colnames(mt))

## Read amplicons 
aa = read.csv(opt$amplicons, h=T, stringsAsFactors=F, sep='\t')

## Select tumor-model pairs in official pairs file 
mt = mt[mt$tumor %in% tn$tumor | mt$model %in% tn$tumor, ]
mt$pair_name = paste0(mt$tumor,'--', mt$model)
mt$in_file = paste0(opt$fc_dir,'/',mt$pair_name,'_feature_similarity_scores.tsv')



## Only take files that exist
message('Subsetting to existing tumor-model comparisons')
mt = mt[file.exists(mt$in_file), ]


## Read input files and merge 
res = do.call(rbind, lapply(mt$in_file, read.csv, h=T, stringsAsFactors=F, sep='\t'))
res$model = gsub('--.*', '', res$Amp1)
res$tumor = gsub('--.*', '', res$Amp2)


## Confirm all Amp1 are tumors and all Amp2 are models
if (!all(res$model %in% mt$model)) {
  stop('Not all Amp1 are model-derived!')
} else {
  model.field = 'Amp1'
}

if (!all(res$tumor %in% mt$tumor)) {
  stop('Not all Amp2 are tumor-derived!')
} else {
  tumor.field = 'Amp2'
}

## For each comparison, pull in tumor and model amplicon-specific input 
col.of.interest = c('Oncogenes', 'Classification', 'Feature.median.copy.number', 'Filter.flag')
res[, paste0('tumor_',tolower(col.of.interest))] = aa[match(res[, tumor.field], aa$Feature.ID), col.of.interest]
res[, paste0('model_',tolower(col.of.interest))] = aa[match(res[, model.field], aa$Feature.ID), col.of.interest]


## Remove anything that was filtered in both the tumor and the model
res = res[res$model_filter.flag == 'None' | res$tumor_filter.flag == 'None', ]


## Make a call on whether these are similar enough to be called concordant 
res$JaccardGenomicSegment_FilterPass = res$JaccardGenomicSegment >= opt$min_jaccard_similarity


## Add unique per-comparison ID to track overlapping amplicons
res$comparison_id = make.unique(paste0(res$tumor,'--',res$model,'_comparison'), sep='_')



##################
## Write result ##
##################

write.table(res, opt$out_file, row.names=F, col.names=T, quote=F, sep='\t')
message(opt$out_file)



##################################
## Filter for concordant ecDNAs ##
##################################

col.sel = c('comparison_id', 'tumor', 'model', 'Amp1', 'Amp2', 'JaccardGenomicSegment', 'AmpOverlapLen', 'Amp1AmpLen', 'Amp2AmpLen', 
            'tumor_oncogenes', 'model_oncogenes', 'tumor_feature.median.copy.number', 'model_feature.median.copy.number', 
            'tumor_filter.flag', 'model_filter.flag', 'JaccardGenomicSegment_FilterPass')

res = res[res$model_classification == 'ecDNA' & res$tumor_classification == 'ecDNA', col.sel]

## Reformat column names 
colnames(res) = gsub('\\.', '_', colnames(res))

## Reformat genes for readability
res$tumor_oncogenes = sapply(parse.genes(res$tumor_oncogenes), paste, collapse=',')
res$model_oncogenes = sapply(parse.genes(res$model_oncogenes), paste, collapse=',')


## Order by overlap 
res = res[order(res$JaccardGenomicSegment, decreasing=T), ]


## Filter for ecDNAs passing similarity filter
res = res[res$JaccardGenomicSegment_FilterPass, ]


## Write table 
out.file.concordant = gsub('\\.txt$', '.concordant_ecdna.txt', opt$out_file)
write.table(res, out.file.concordant, row.names=F, col.names=T, quote=F, sep='\t')
message(out.file.concordant)
