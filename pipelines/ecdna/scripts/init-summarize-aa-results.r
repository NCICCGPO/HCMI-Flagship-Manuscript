#!/nfs/sw/R/R-4.0.0/bin/Rscript
## Summarize AmpliconArchitect results
libs = c('optparse', 'ggplot2', 'ggbeeswarm')
invisible(suppressPackageStartupMessages(sapply(libs, require, character.only=T)))
options(width=200, scipen=999)



## Read finish flag, handling missing files
read.finish.flag = function(f) {

  sapply(f, function(x) ifelse(file.exists(x), 
                               readLines(x)[1], 
                               'UNSUCCESSFUL'))

}



## Get arguments
option_list = list(
  make_option(c("-t", "--tn_file"),                          type='character', help="Tumor-normal pairs file"),
  make_option(c("-T", "--tm_file"),                          type='character', help="Tumor-model pairs file"),
  make_option(c("-a", "--aa_dir"),                           type='character', help="Top-level AmpliconArchitect output directory"),
  make_option(c("-o", "--out_file"),                         type='character', help="Output PDF"),
  make_option(c("-O", "--out_file_txt"),                     type='character', help="Output table (filtered)"),
  make_option(c("-u", "--out_file_txt_unfiltered"),          type='character', help="Output table (UNfiltered)"))
opt = parse_args(OptionParser(option_list=option_list))



## Read tumor normal pairs 
tn = read.csv(opt$tn_file, h=T, stringsAsFactors=F)
tn$tumor = tn$tumor_aliquot_barcode
tn$normal = tn$normal_aliquot_barcode


# Determine whether runs were successful or not 
## Only consider TN pairs with successful runs 
tn$finish_flag = paste0(opt$aa_dir,'/', tn$participantId, '/', tn$pairId, '/', tn$pairId,'_finish_flag.txt')
tn$finish_flag = read.finish.flag(tn$finish_flag)


tn = tn[tn$finish_flag == 'All stages completed', ]

dim(tn)


## Read tumor-model pairs 
mt = read.csv(opt$tm_file, h=T, stringsAsFactors=F, sep='\t')
colnames(mt) = tolower(colnames(mt))
tn$model = tn$tumor %in% mt$model

n.model = sum(tn$model)
n.tumor = sum(!tn$model)


## Read AA amplicon classification summary
tn$amplicon_classification_dir = paste0(opt$aa_dir,'/',tn$participantId,'/',tn$pairId,'/',tn$pairId,'_classification')
tn$amplicon_classification_profiles = paste0(tn$amplicon_classification_dir,'/',tn$pairId,'_result_table.tsv')

head(tn$amplicon_classification_profiles)

missing.profiles = which(!file.exists(tn$amplicon_classification_profiles))
if (length(missing.profiles) > 0) {
  message('Missing amplicon classification profiles for ', length(missing.profiles),'/',nrow(tn),' samples. Skipping...')
  tn = tn[-missing.profiles, ]
}

aa = lapply(tn$amplicon_classification_profiles, read.csv, h=T, stringsAsFactors=F, sep='\t')
aa = Reduce(rbind, aa)
aa$Filter.flag[is.na(aa$Filter.flag)] = 'None'

## Only filter for amplicons with a classification
aa = aa[!is.na(aa$Classification), ]
aa$model = tn$model[match(aa$Sample.name, tn$pairId)]


## Remove anything with a filter
aa.unfiltered = aa
aa = aa[aa$Filter.flag == 'None', ]


## Overall count summary
aa.event.counts = as.data.frame(table(aa$Classification))
colnames(aa.event.counts) = c('event_type', 'count')
aa.event.counts



###################
## Tumor summary ##
###################

message('\n-- Tumor Summary --')

num.samples.with.amplicon = length(unique(aa$Sample.name[!aa$model]))
pct.samples.with.amplicon = round(num.samples.with.amplicon / n.tumor, 4) * 100
message('Samples with an amplicon detected: ',num.samples.with.amplicon,'/',n.tumor, ' (', pct.samples.with.amplicon,'%)')

num.samples.with.cyclic = length(unique(aa$Sample.name[aa$Classification == 'ecDNA' & !aa$model]))
pct.samples.with.cyclic = round(num.samples.with.cyclic / n.tumor, 4) * 100
message('Samples with a ecDNA amplicon detected: ',num.samples.with.cyclic,'/',n.tumor, ' (', pct.samples.with.cyclic,'%)')

num.samples.with.cpx.noncyclic = length(unique(aa$Sample.name[aa$Classification == 'Complex-non-cyclic' & !aa$model]))
pct.samples.with.cpx.noncyclic = round(num.samples.with.cpx.noncyclic / n.tumor, 4) * 100
message('Samples with a COMPLEX NON-CYCLIC amplicon detected: ',num.samples.with.cpx.noncyclic,'/',n.tumor, ' (', pct.samples.with.cpx.noncyclic,'%)')

num.samples.with.linear.amp = length(unique(aa$Sample.name[aa$Classification == 'Linear' & !aa$model]))
pct.samples.with.linear.amp = round(num.samples.with.linear.amp / n.tumor, 4) * 100
message('Samples with a LINEAR AMPLIFICATION amplicon detected: ',num.samples.with.linear.amp,'/',n.tumor, ' (', pct.samples.with.linear.amp,'%)')

message('\n\n')


message('\n-- Model Summary --')

num.samples.with.amplicon = length(unique(aa$Sample.name[aa$model]))
pct.samples.with.amplicon = round(num.samples.with.amplicon / n.model, 4) * 100
message('Samples with an amplicon detected: ',num.samples.with.amplicon,'/',n.model, ' (', pct.samples.with.amplicon,'%)')

num.samples.with.cyclic = length(unique(aa$Sample.name[aa$Classification == 'ecDNA' & aa$model]))
pct.samples.with.cyclic = round(num.samples.with.cyclic / n.model, 4) * 100
message('Samples with a ecDNA amplicon detected: ',num.samples.with.cyclic,'/',n.model, ' (', pct.samples.with.cyclic,'%)')

num.samples.with.cpx.noncyclic = length(unique(aa$Sample.name[aa$Classification == 'Complex-non-cyclic' & aa$model]))
pct.samples.with.cpx.noncyclic = round(num.samples.with.cpx.noncyclic / n.model, 4) * 100
message('Samples with a COMPLEX NON-CYCLIC amplicon detected: ',num.samples.with.cpx.noncyclic,'/',n.model, ' (', pct.samples.with.cpx.noncyclic,'%)')

num.samples.with.linear.amp = length(unique(aa$Sample.name[aa$Classification == 'Linear' & aa$model]))
pct.samples.with.linear.amp = round(num.samples.with.linear.amp / n.model, 4) * 100
message('Samples with a LINEAR AMPLIFICATION amplicon detected: ',num.samples.with.linear.amp,'/',n.model, ' (', pct.samples.with.linear.amp,'%)')

message('\n\n')




## Summarize per-sample counts
# message('Amplicon per-sample count summary')
# summary(as.numeric(table(aa$pairId[aa$amplicon_decomposition_class != 'No amp/Invalid'])))

message('ecDNA per-sample count summary')
summary(as.numeric(table(aa$Sample.name[aa$Classification == 'ecDNA'])))
tapply(aa$Sample.name[aa$Classification == 'ecDNA'], 
       aa$model[aa$Classification == 'ecDNA'], function(x) summary(as.numeric(table(x))))

# sort(table(aa$Sample.name[aa$Classification == 'ecDNA']))

# message('COMPLEX NON-CYCLIC per-sample count summary')
# summary(as.numeric(table(aa$Sample.name[aa$Classification == 'Complex-non-cyclic'])))

# message('LINEAR AMPLIFICATION per-sample count summary')
# summary(as.numeric(table(aa$Sample.name[aa$Classification == 'Linear'])))



# ## Confirm that all Cyclic calls are ecDNA+ 
# if (!all(aa$`ecDNA+`[aa$amplicon_decomposition_class == 'Cyclic'] == 'Positive')) {
#   stop('ERROR: Not all cyclic amplicons are ecDNA+!')
# } else {
#   message('All cyclic amplicons are ecDNA+!')
# }


# ## Check to see if we detected any BFBs
# if (any(aa$`BFB+` == 'Positive')) {
#   stop('AmpliconArchitect detected BFB cycles!')
# } else {
#   message('No BFB cycles detected!')
# }






## Size summary
aa$Captured.interval.length.mb = aa$Captured.interval.length / 1E6

table(aa$Classification)

message('Models with ecDNA: ', length(unique(aa$Sample.name[aa$model & aa$Classification == 'ecDNA'])))
message('Tumors with ecDNA: ', length(unique(aa$Sample.name[!aa$model & aa$Classification == 'ecDNA'])))


aa$Classification = factor(aa$Classification, names(sort(tapply(aa$Captured.interval.length.mb, aa$Classification, median, na.rm=T))))


message('\n-- Size summary --')
med.ecdna.size = round(median(aa$Captured.interval.length.mb[aa$Classification == 'ecDNA']), 3)
min.ecdna.size = round(min(aa$Captured.interval.length.mb[aa$Classification == 'ecDNA']), 3)
max.ecdna.size = round(max(aa$Captured.interval.length.mb[aa$Classification == 'ecDNA']), 3)

message('Median size across cohort: ',med.ecdna.size,'MB (', min.ecdna.size ,'-',max.ecdna.size,')')

message('Size summary by model (TRUE) and tumor (FALSE)')
tapply(aa$Captured.interval.length.mb[aa$Classification == 'ecDNA'], aa$model[aa$Classification == 'ecDNA'], summary)

message('Size comparison:')
wilcox.test(Captured.interval.length.mb ~ model, data=aa[aa$Classification == 'ecDNA', ])



## CN summary
message('\n-- CN summary --')
med.ecdna.cn = round(median(aa$Feature.median.copy.number[aa$Classification == 'ecDNA']), 3)
min.ecdna.cn = round(min(aa$Feature.median.copy.number[aa$Classification == 'ecDNA']), 3)
max.ecdna.cn = round(max(aa$Feature.median.copy.number[aa$Classification == 'ecDNA']), 3)

message('Median CN across cohort: ',med.ecdna.cn,' (', min.ecdna.cn ,'-',max.ecdna.cn,')')

message('CN summary by model (TRUE) and tumor (FALSE)')
tapply(aa$Feature.median.copy.number[aa$Classification == 'ecDNA'], aa$model[aa$Classification == 'ecDNA'], summary)

message('CN comparison:')
wilcox.test(Feature.median.copy.number ~ model, data=aa[aa$Classification == 'ecDNA', ])



## Gene summary
message('\n-- Oncogene summary --')
oncogenes = gsub('\\]', ')', gsub('\\[', 'c(', aa$Oncogenes[aa$Classification == 'ecDNA']))
oncogenes = lapply(oncogenes, function(x) eval(parse(text=x)))
write.table(head(as.data.frame(sort(table(unlist(oncogenes)), decreasing=T)), 20), row.names=F, quote=F)


message('\n-- Oncogene summary : TUMORS --')
oncogenes = gsub('\\]', ')', gsub('\\[', 'c(', aa$Oncogenes[aa$Classification == 'ecDNA' & !aa$model]))
oncogenes = lapply(oncogenes, function(x) eval(parse(text=x)))
write.table(head(as.data.frame(sort(table(unlist(oncogenes)), decreasing=T)), 20), row.names=F, quote=F)


message('\n-- Oncogene summary : MODELS --')
oncogenes = gsub('\\]', ')', gsub('\\[', 'c(', aa$Oncogenes[aa$Classification == 'ecDNA' & aa$model]))
oncogenes = lapply(oncogenes, function(x) eval(parse(text=x)))
write.table(head(as.data.frame(sort(table(unlist(oncogenes)), decreasing=T)), 20), row.names=F, quote=F)


###########################
## Write merged summmary ##
###########################

write.table(aa.unfiltered, opt$out_file_txt_unfiltered, row.names=F, col.names=T, quote=F, sep='\t')
message(opt$out_file_txt_unfiltered)

write.table(aa, opt$out_file_txt, row.names=F, col.names=T, quote=F, sep='\t')
message(opt$out_file_txt)



##############
## Plotting ##
##############

# pdf(opt$out_file, width=5, height=5)
pdf(opt$out_file, width=6, height=4)

## Event size 
aa.class.cols = c(`ecDNA`='#E07694', `Linear`='#2C8EC1', `BFB`='firebrick', `Complex`='#45C166')
levels(aa$Classification)[levels(aa$Classification) == 'Complex-non-cyclic'] = 'Complex'

ggplot(aa, aes(x=Classification, y=Captured.interval.length.mb , fill=Classification)) +
  geom_violin() +
  # geom_quasirandom(alpha=0.5) +
  geom_boxplot(width=0.15, outlier.shape=NA, fill='#FFFFFF') +
  scale_fill_manual(values=aa.class.cols) +
  scale_y_log10() + 
  xlab('AmpliconArchitect classification') +
  ylab('Event size (MB)') +
  theme_bw() +
  theme(text=element_text(size=19),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank(),
        legend.position='none')



## ecDNAs per sample
ecdnas.per.sample = tapply(aa$Sample.name[aa$Classification == 'ecDNA'],
                      aa$model[aa$Classification == 'ecDNA'], 
                      function(x) unname(table(x)))

## Check that ecDNAs per sample aren't different between groups 
wilcox.test(ecdnas.per.sample[[1]], ecdnas.per.sample[[2]])

ecdnas.per.sample = data.frame(`ecDNA count`=unlist(ecdnas.per.sample), check.names=F)
count.states = sort(unique(ecdnas.per.sample$`ecDNA count`))

ggplot(ecdnas.per.sample, aes(x=`ecDNA count`)) +
  geom_histogram(binwidth=1, color='#FFFFFF') + 
  scale_x_continuous(name='ecDNAs per sample', breaks=count.states, labels=count.states) +
  ylab('Sample count') +
  theme_bw() +
  theme(text=element_text(size=17),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank()) 



## Median event CN
aa$Classification = factor(aa$Classification, names(sort(tapply(aa$Feature.median.copy.number, aa$Classification, median, na.rm=T))))

ggplot(aa, aes(x=Classification, y=Feature.median.copy.number, fill=Classification)) +
  geom_violin() +
  # geom_quasirandom() +
  geom_boxplot(width=0.15, outlier.shape=NA, fill='#FFFFFF') +
  scale_fill_manual(values=aa.class.cols) +
  scale_y_log10() + 
  xlab('AmpliconArchitect classification') +
  ylab('Median CN') +
  theme_bw() +
  theme(text=element_text(size=19),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank(),
        legend.position='none')




dev.off()
message(opt$out_file)
