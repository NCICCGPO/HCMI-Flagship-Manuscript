## Build input tables for joint calling 
libs = c('optparse')
invisible(suppressPackageStartupMessages(sapply(libs, require, character.only=T)))
options(width=200, scipen=999)


## Get arguments
option_list = list(
  make_option(c("-t", "--tn_file"),  type='character', help="Tumor-normal pairs CSV"),
  make_option(c("-c", "--cnv_dir"),  type='character', help="Dir with CNV profiles"),
  make_option(c("-b", "--bam_map"),  type='character', help="Map between aliquot IDs and BAM paths"),
  make_option(c("-o", "--out_dir"),  type='character', help="Output directory"))
opt = parse_args(OptionParser(option_list=option_list))


## Read tumor normal pairs and bam map 
tn = read.csv(opt$tn_file, h=T, stringsAsFactors=F)
bam = read.csv(opt$bam_map, h=T, stringsAsFactors=F)

## Map BAMs to aliquots 
tn$bam = bam$aliquot_barcode_named_path[match(tn$tumor_aliquot_barcode, bam$aliquot_barcode)]

## Generate CNV paths 
tn$cnv = paste0(opt$cnv_dir, '/', tn$pairId,'.jabba.bed')

## Filter for samples that we have CNV profiles for 
tn = tn[file.exists(tn$cnv), ]

## Reformat to what AA expects 
tn$specimen = 'tumor'
tn = tn[, c('participantId', 'pairId', 'bam', 'specimen', 'cnv')]

## Split on participant
tn = split(tn, tn$participantId)



## Write out one input file per participant 
for (i in 1:length(tn)) {

  out.file = paste0(opt$out_dir,'/', names(tn)[i],'.joint_aa_input.txt')
  write.table(tn[[i]][, -1], out.file, row.names=F, col.names=F, quote=F, sep='\t')
  message(out.file)

}
