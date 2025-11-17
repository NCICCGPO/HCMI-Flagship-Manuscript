#!/usr/bin/env python

# Script to generate bamreadcount like features from tumor and normal bam files
####################################
# Version 1 (2017-06-22)
# Kanika Arora (karora@nygenome.org)
# Rashesh Sanghvi (rsanghvi@nygenome.org)
# New York Genome Center
####################################

import pysam
import sys
import argparse
import os
import subprocess

class SomaticSNVIndel:
        
    def __init__(self, chr, pos, id, ref, alt):
        ''' Object to hold single variant'''
        self.chr = chr
        self.pos = pos
        self.id = id
        self.ref = ref
        self.alt = alt
        self.tumor = AlleleCountFeatures(prefix="T")
        self.normal = AlleleCountFeatures(prefix="N")
    
    @staticmethod
    def print_header():
        tumor = AlleleCountFeatures(prefix="T")
        normal = AlleleCountFeatures(prefix="N")
        header = ['#Chr','Pos','ID','Ref','Alt']
        header += tumor.get_header()
        header += normal.get_header()
        return header
    
    def print_variant(self):
        arr=[self.chr, self.pos, self.id, self.ref, self.alt]
        arr+=self.tumor.get_bamfeatures()
        arr+=self.normal.get_bamfeatures()
        return map(str,arr)


class AlleleCountFeatures:
    def __init__(self, prefix):
        ''' Object to hold single allele within a position'''
        self.prefix=prefix
        self.dp = 0
        self.uniq_pairs_dp = 0
        self.uniq_alt_cnt = 0
        self.uniq_ref_cnt = 0
        self.uniq_other_cnt = 0
        self.ref_cnt = 0
        self.alt_cnt = 0
        self.other_cnt = 0
        self.allels={}

    def get_header(self):
        prefix=self.prefix
        header=['dp','uniq_pairs_dp',
                'ref_cnt','alt_cnt',
                'VAF','other_cnt', 
                'uniq_ref_cnt', 'uniq_alt_cnt', 
                'uniq_VAF', 'uniq_other_cnt', 
                'Alleles']
        header = [prefix + "_" + ele for ele in header]
        return header

    def get_bamfeatures(self):
        '''Use counts from BAM to calculate allele counts'''
        ##nor_ref, nor_alt, tum_ref, tum_alt
        features=[]
        if self.dp == 0:
            VAF = 0.0
        else:
            VAF = float(self.alt_cnt)/float(self.dp)
        if self.uniq_pairs_dp == 0:
            uniq_VAF = 0.0
        else:
            uniq_VAF = self.uniq_alt_cnt/float(self.uniq_pairs_dp)
        features += [self.dp, self.uniq_pairs_dp, 
                     self.ref_cnt, self.alt_cnt, 
                     VAF, self.other_cnt, self.uniq_ref_cnt, 
                     self.uniq_alt_cnt, uniq_VAF, 
                     self.uniq_other_cnt]
        tmp_ft = ""
        for alleles in self.allels.keys():
            tmp_ft += alleles+":" + str(self.allels[alleles]) + ";"
        features += [tmp_ft]
        return features

##########################################################################################
##########################################################################################

class ArgumentParser(argparse.ArgumentParser):
    def error(self, message):
        self.print_help(sys.stderr)
        self.exit(2, '\nERROR: %s\n\n' % (message))

def get_args():
    parser = ArgumentParser(prog='parse_bam_generate_features', description='Parses bam file and generates features for positions in VCF file', epilog='', formatter_class=lambda prog: argparse.ArgumentDefaultsHelpFormatter(prog, max_help_position=100, width=150))
    parser.add_argument('--tumor_bam', help='Tumor bam file to be parsed.', required=True)
    parser.add_argument('--normal_bam', help='Normal bam file to be parsed.', required=True)
    parser.add_argument('--vcf', help = 'SNV VCF file.', required=True)
    parser.add_argument('--chrom', help = 'Chrom to resrict run to.', default=False)
    parser.add_argument('--output', help = 'Output txt file.', required=True)
    parser.add_argument('--reference', help='Reference genome FASTA file', required=True)
    parser.add_argument('-i', '--max_indel_len_for_count',
                            help='Maximum indel length for generating counts',
                            default=10, type=int)

    parser.add_argument('--min_base_quality', help='Minimum base quality', default=10, type=int)
    parser.add_argument('--min_mapping_quality', help='Minimum mapping quality', default=10, type=int)
    parser.add_argument('--no_header', help='Do not add header line to output.', action='store_true')
    args=parser.parse_args()
    return args


def get_type(ref, alt):
    '''
        Fill in the type field.
    '''
    if len(ref) == 1 and len(alt) == 1:
        type = 'SNV'
    elif ',' in alt:
        type = 'MULTI'
    elif len(ref) == 1 and len(alt) > 1 and ref[0] == alt[0]:
        type = 'INS'
    elif len(ref) > 1 and len(alt) == 1 and ref[0] == alt[0]:
        type = 'DEL'
    elif len(ref) > 1 \
        and len(ref) == len(alt):
        type = 'MNV'
    else:
        type = 'COMPLEX'
    return type


def is_complex(ref, alt):
    '''
        Test if a variant is complex. Skips variants
        with an anchor base.
    '''
    complex = False
    if len(ref) != len(alt) \
        and len(ref) > 1 \
        and len(alt) > 1 \
        and ref[0] != alt[0]:
            complex = True
    return complex


def is_too_long(ref, alt, MAX_INDEL_LEN):
    '''
    Test if an INDEL is too long given the length cut off.
    '''
    too_long = False
    if (len(ref) > MAX_INDEL_LEN or len(alt) > MAX_INDEL_LEN):
        too_long = True
    return too_long


def parse_pileup(pileup, SNP, chr, pos, ref, alt, ReferenceFasta,
                 MIN_BQ=10, MIN_MQ=10, MAX_INDEL_LEN=10):
    '''Translate piplup into allele counts'''
    alt_read_names = dict()
    ref_read_names = dict()
    other_read_names = dict()
    #### Skip unsupported variants#######
    too_long = is_too_long(ref, alt, MAX_INDEL_LEN)
    complex = is_complex(ref, alt)
    if complex or too_long:
        print("%i\t%s\t%s  This type of variants are currently Not supported (long, complex)." %(pos, ref, alt))
        return
    kind = get_type(ref, alt)
    if kind in ['MULTI', 'MNV', 'COMPLEX']:
        print("%i\t%s\t%s  This type of variants are currently Not supported." %(pos, ref, alt))
        return

    og_alt=alt
    og_ref=ref
    #### Incorporating INDEL calling -rsanghvi#######
    for pileupcolumn in pileup:
        if pileupcolumn.pos != pos:
            continue
        for pileupread in pileupcolumn.pileups:
            # if ref position is N take the next position
            if pileupread.is_refskip:
                continue
            # if the position in the read is .is_del pos is none so take next
            pos_in_read = pileupread.query_position
            # skip reads where the position is already a deletion (is_del)
            if not pos_in_read:
                continue
            BQ = pileupread.alignment.query_qualities[pos_in_read]
            MQ = pileupread.alignment.mapping_quality
            if BQ < MIN_BQ \
                    or MQ < MIN_MQ \
                    or pileupread.alignment.is_supplementary:
                continue
            base = pileupread.alignment.query_sequence[pos_in_read]
            alt_allele = base
            indel = pileupread.indel
            align_start = pileupread.alignment.query_alignment_start
            align_end = pileupread.alignment.query_alignment_end
            align_len = pileupread.alignment.query_alignment_length
            proper_paired = pileupread.alignment.is_proper_pair
            read_name = pileupread.alignment.query_name
            read_length = pileupread.alignment.query_length
            pos_as_frac_read = min(read_length - pos_in_read - 1, pos_in_read)/(float(read_length)/2) ## if base is at center of read, this number will be 1, if base is at end of read, this number will be 0
            pos_as_frac_align = min(pos_in_read - align_start, align_end - pos_in_read - 1)/(float(align_len)/2) ## if base is at center of aligned part of the read, this number will be 1, if base is at end of alignment, this number will be 0
            cigtup = pileupread.alignment.cigartuples
            cig = pileupread.alignment.cigarstring
            NM = pileupread.alignment.get_tag('NM') ## gives edit distance.
            frac_edit_dist = NM/align_len
            SNP.dp+=1

            #### Incorporating INDEL calling -rsanghvi#######			
            if (kind == 'DEL'):     
                #### TO MATCH THE NEW BAMREADCOUNT. The counts are being made from the first base in deletion and not from the Anchor position.
                ref = og_ref[:2]
                alt = "-"+og_ref
                base = pileupread.alignment.query_sequence[pos_in_read:pos_in_read+2]
                alt_allele = base
            if (kind == 'INS'):
                ref = og_ref
                alt = "+"+og_alt	
            ##Change the base sequence if there is an INDEL in the read. And count appropriately depending on what the force call mutation is.
            if indel < 0:
                base = "-"+pileupread.alignment.query_sequence[pos_in_read] + ReferenceFasta.fetch(reference=chr,
                                                                                                   start=pos+1,
                                                                                                   end=pos+abs(indel)+1) ### This represents the sequence that is deleted
                alt_allele="-"+ReferenceFasta.fetch(reference=chr, start=pos+1, end=pos+abs(indel)+1)
                if (kind != 'DEL'): 
                    ## If the force call mutation is not deletion, but the read has a deletion, we currently record the deletion in Allele information, but just check the position at hand for ref or alt base
                    base = pileupread.alignment.query_sequence[pos_in_read]
                    if pileupread.alignment.query_sequence[pos_in_read] in SNP.allels:
                        SNP.allels[pileupread.alignment.query_sequence[pos_in_read]]+=1
                    else:
                        SNP.allels[pileupread.alignment.query_sequence[pos_in_read]]=1
            elif indel > 0: 
                base = "+"+pileupread.alignment.query_sequence[pos_in_read:pos_in_read+abs(indel)+1]   
                ### This represents the sequence after the exact insertion in the White list.
                alt_allele="+"+pileupread.alignment.query_sequence[pos_in_read+1:pos_in_read+abs(indel)+1]
                if (kind != 'INS'):  ## If the force call mutation is not insertion, but the read has an insertion, we currently record the deletion in Allele information, but just check the position at hand for ref or alt base
                    base = pileupread.alignment.query_sequence[pos_in_read]
                    if pileupread.alignment.query_sequence[pos_in_read] in SNP.allels:
                        SNP.allels[pileupread.alignment.query_sequence[pos_in_read]]+=1
                    else:
                        SNP.allels[pileupread.alignment.query_sequence[pos_in_read]]=1
            #### Incorporating INDEL calling -rsanghvi#######
            if alt_allele in SNP.allels:
                SNP.allels[alt_allele]+=1
            else:
                SNP.allels[alt_allele]=1
            if base == ref:
                ref_read_names[read_name]=1
                SNP.ref_cnt+=1
            elif base == alt:
                alt_read_names[read_name]=1
                SNP.alt_cnt+=1
            else:
                other_read_names[read_name]=1
                SNP.other_cnt+=1
    for key in list(ref_read_names.keys()):
        if key in alt_read_names:
            del alt_read_names[key]
            del ref_read_names[key]
        if key in other_read_names:
            del other_read_names[key]
            del ref_read_names[key]
    for key in list(alt_read_names.keys()):
        if key in other_read_names:
            del alt_read_names[key]
            del other_read_names[key]

    SNP.uniq_pairs_dp = len(alt_read_names.keys()) + len(ref_read_names.keys()) + len(other_read_names.keys())
    SNP.uniq_ref_cnt = len(ref_read_names.keys())
    SNP.uniq_alt_cnt = len(alt_read_names.keys())
    SNP.uniq_other_cnt = len(other_read_names.keys())
    return SNP

# Load bam files
def reset_token():
        try:
            result = subprocess.run(['gcloud','auth', 'application-default', 'print-access-token'],
                                    check=True,
                                    stdout=subprocess.PIPE).stdout.decode('utf-8')
            os.environ['GCS_OAUTH_TOKEN'] = result.rstrip()
        except subprocess.CalledProcessError as err:
            log.error(err.output.decode('utf-8'))
            log.error('Failed to regenerate an access token')
            return False
        return True


def open_bam(BAM_FILE, initial=True):
    try:
        samfile = pysam.AlignmentFile(BAM_FILE, "rb")
    except (PermissionError, AssertionError):
        if initial:
            reset_token()
            samfile = open_bam(BAM_FILE, initial=False)
        else:
            log.error('PermissionError and print-access-token failed to recover')
            sys.exit(1)
    return samfile


def fetch_pileups(chrom, pos, samfile, BAM_FILE, initial=True):
    try:
        pileup = samfile.pileup(chrom, int(pos)-1, int(pos),
                                max_depth=100000)
        return pileup
    except (PermissionError, AssertionError):
        if initial:
            reset_token()
            samfile = open_bam(BAM_FILE, initial=True)
            pileup = samfile.pileup(chrom, int(pos)-1, int(pos),
                        max_depth=100000)
        else:
            log.error('PermissionError and print-access-token failed to recover')
            sys.exit(1)


def process_variants(args, OUT, VCF, chrom, TBAM, NBAM,
                    TumorSamFile, NormalSamFile,
                    MAX_INDEL_LEN, ReferenceFasta):
    ''' Process variants '''
    seen_header = 0
    if args.no_header is True:
        seen_header = 1
    with open(OUT, "w") as o:
        # print header (even if VCF has no variants)
        o.write("\t".join(SomaticSNVIndel.print_header())  + "\n")
        with open(VCF) as f:
            for line in f:
                line = line.strip()
                # skip comments
                if line.startswith("#"):
                    continue
                 # write out header
                toks = line.split()
                # skip non-chrom
                if chrom:
                    if toks[0] != chrom:
                        continue
                # begin storing variants for printing
                if seen_header == 0:
                    seen_header=1
                # load variant object
                variant = SomaticSNVIndel(chr=toks[0], pos=toks[1],
                                          id=toks[2], ref=toks[3], alt=toks[4])
                # generate pileups
                T_pileup = fetch_pileups(chrom=toks[0],
                                        pos=int(toks[1]),
                                        samfile=TumorSamFile,
                                        BAM_FILE=TBAM,
                                        initial=True)
                variant.tumor = parse_pileup(T_pileup, variant.tumor,
                                             toks[0], int(toks[1])-1, toks[3], toks[4],
                                             ReferenceFasta=ReferenceFasta,
                                             MIN_BQ=args.min_base_quality,
                                             MIN_MQ=args.min_mapping_quality,
                                             MAX_INDEL_LEN=MAX_INDEL_LEN)
                N_pileup = fetch_pileups(chrom=toks[0],
                                    pos=int(toks[1]),
                                    samfile=NormalSamFile,
                                    BAM_FILE=NBAM,
                                    initial=True)
                variant.normal = parse_pileup(N_pileup, variant.normal,
                                            toks[0], int(toks[1])-1, toks[3], toks[4],
                                            ReferenceFasta=ReferenceFasta,
                                            MIN_BQ=args.min_base_quality,
                                            MIN_MQ=args.min_mapping_quality,
                                            MAX_INDEL_LEN=MAX_INDEL_LEN)
                # write out
                o.write("\t".join(variant.print_variant()) + "\n")


def run(TBAM, NBAM, VCF,
        REF, chrom, args,
        MAX_INDEL_LEN, OUT):
    TumorSamFile = open_bam(TBAM, initial=True)
    NormalSamFile = open_bam(NBAM, initial=True)
    ReferenceFasta = pysam.FastaFile(REF)
    process_variants(args, OUT, VCF, chrom, TBAM, NBAM,
                    TumorSamFile, NormalSamFile,
                    MAX_INDEL_LEN, ReferenceFasta)


def main():
    args = get_args()
    MAX_INDEL_LEN=args.max_indel_len_for_count
    TBAM=args.tumor_bam
    NBAM=args.normal_bam
    REF=args.reference
    VCF=args.vcf
    OUT=args.output
    run(TBAM, NBAM, VCF,
        REF, args.chrom, args,
        MAX_INDEL_LEN, OUT)


if __name__ == "__main__":
    main()
