## Extract CNV calls from JaBbA VCF to BED, retaining only copy number info 
libs = c('optparse', 'reshape2', 'stringr', 'gUtils')
invisible(suppressPackageStartupMessages(sapply(libs, require, character.only=T)))
options(width=200)


VCF_COLS     = c('CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'DATA')
BED_COLS     = c('CHROM', 'START', 'END', 'CN')
BED_COLS_OUT = c('#chr', 'start', 'end', 'cn')
GR_COLS      = c('seqnames', 'start', 'end', 'CN', 'gene')
GR_COLS_OUT  = c('#chr', 'start', 'end', 'cn', 'gene')


## Get arguments
option_list = list(
  make_option(c("-i", "--in_file"),   type='character', help="Input VCF (jabba.simple.cnv.vcf)"),
  make_option(c("-g", "--gene_bed"),  type='character', help="Optional BED file with gene coordinates"),
  make_option(c("-o", "--out_file"),  type='character', help="Output BED"))
opt = parse_args(OptionParser(option_list=option_list))


## Read VCF
vcf = read.csv(opt$in_file, h=F, stringsAsFactors=F, sep='\t', comment.char='#')
colnames(vcf) = VCF_COLS


## Extract start, end, copy number
vcf$START = str_extract(vcf$INFO, '(?<=START=)[0-9]+(?=;)')
vcf$END = str_extract(vcf$INFO, '(?<=END=)[0-9]+(?=;)')
vcf$CN = str_extract(vcf$DATA, '(?<=:)[0-9]+')

vcf = vcf[, BED_COLS]


## If we have a gene BED file 
if (!is.null(opt$gene_bed)) {
  
  ## Read genes
  genes = read.csv(opt$gene_bed, h=F, stringsAsFactors=F, sep='\t')
  colnames(genes) = c('chr','start','end','gene')
  genes = makeGRangesFromDataFrame(genes, keep.extra.columns=T)
  
  ## Aggregate genes across ranges in vcf 
  vcf = makeGRangesFromDataFrame(vcf, keep.extra.columns=T)
  vcf = as.data.frame(vcf %$% genes)
  vcf$gene = gsub(', ', ',', vcf$gene, fixed=T)
  
  vcf = vcf[, GR_COLS]
  colnames(vcf) = GR_COLS_OUT
  
} else {
  
  colnames(vcf) = BED_COLS_OUT
  
}


## Write to disk
write.table(vcf, opt$out_file, col.names=T, row.names=F, sep='\t', quote=F)
