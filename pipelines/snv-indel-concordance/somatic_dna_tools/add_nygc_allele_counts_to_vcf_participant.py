#!/usr/bin/env python
#	USAGE: python add_nygc_allele_counts_to_vcf_participant.py --help
#   DESCRIPTION: Runs pileup on tumor(s) and normal bam files to compute allele counts for 
#   bi-allelic SNV and Indel variants in VCF file and adds pileup format entries to the VCF file.
################################################################################
##################### COPYRIGHT ################################################
# New York Genome Center

# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright (2024) by the New York
# Genome Center. All rights are reserved. This software is supplied without
# any warranty or guaranteed support whatsoever. The New York Genome Center
# cannot be responsible for its use, misuse, or functionality.
# Version: 0.1
# Author: Kanika Arora (karora@nygenome.org) and Jennifer M Shelton
##################### /COPYRIGHT ###############################################
################################################################################
import argparse
import pysam


def modify_header(bcf_in):
    '''
        Add new FORMAT fields
    '''
    # force count reads
    formats = [
            'aliquot_AD', 
            'aliquot_DP', 
            'aliquot_AF', 
            'BQ',
            'MQ'
        ]
    # ##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles computed using a custom pileup-based approach.">
    bcf_in.header.formats.add(id='AD', 
                           number='R',
                           type='Integer',
                           description='Allelic depths for the ref and alt alleles computed using a custom pileup-based approach (Minimum mapping quality = 10, Minimum base quality = 10).')
    # ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth of coverage: Number of reads covering site computed using a custom pileup-based approach.">
    bcf_in.header.formats.add(id='DP', 
                           number='1',
                           type='Integer',
                           description='Depth of coverage: Number of reads covering site computed using a custom pileup-based approach (Minimum mapping quality = 10, Minimum base quality = 10).')
    # ##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele fraction of alternate allele in the sample computed using a custom pileup-based approach.">
    bcf_in.header.formats.add(id='AF', 
                           number='A',
                           type='Float',
                           description='Allele fraction of alternate allele in the sample computed using a custom pileup-based approach (Minimum mapping quality = 10, Minimum base quality = 10).')
    bcf_in.header.formats.add(id='NORMAL_AF', 
                           number='1',
                           type='Float',
                           description='Allele fraction of reference allele in the sample computed using a custom pileup-based approach (Minimum mapping quality = 10, Minimum base quality = 10).')
    bcf_in.header.formats.add(id='BQ', 
                       number='.',
                       type='Integer',
                       description='Base quality for all reads supporting the Alt at the start position. No min MapQ or Base Quality')
    bcf_in.header.formats.add(id='MQ', 
                       number='.',
                       type='Integer',
                       description='Map quality for all reads supporting the Alt at the start position. No min MapQ or Base Quality')
    return bcf_in

def read_vcf(vcf_file):
    '''
        Read in annotated VCF file.
    '''
    bcf_in = pysam.VariantFile(vcf_file)  # auto-detect input format
    return bcf_in

def read_bam(bam_file):
    bam = pysam.AlignmentFile(bam_file, "rb")
    return bam

def get_bam_sample(bam):
    samples = [bam.header.as_dict()['RG'][i]['SM'] for i in range(len(bam.header.as_dict()['RG']))]
    assert len(set(samples)) == 1, 'BAM must come from one aliquot'
    return list(set(samples))[0]

def is_too_long(ref, alt, MAX_INDEL_LEN):
    '''
    Test if an INDEL or COMPLEX event is too long for computing allele counts using NYGC's pileup method given the length cut off.
    '''
    too_long = False
    if (len(ref) > MAX_INDEL_LEN or len(alt) > MAX_INDEL_LEN):
        too_long = True
    return too_long

def parse_snv_mnv(read_name, pileupread, ref, alt, pos_in_read, anchor_mismatch,
                ref_reads, alt_reads, other_reads, ):
    is_alt = False
    if pileupread.alignment.query_sequence[pos_in_read:pos_in_read + len(ref)] == ref:
        ref_reads.append(read_name)
    elif pileupread.alignment.query_sequence[pos_in_read:pos_in_read + len(alt)] == alt:
        alt_reads.append(read_name)
        if pileupread.indel != 0:
            anchor_mismatch+=1
        is_alt = True
    else:
        other_reads.append(read_name)
    return ref_reads, alt_reads, other_reads, anchor_mismatch, is_alt


def parse_indel(read_name, pileupread, ref, alt, pos_in_read, anchor_mismatch,
                ref_reads, alt_reads, other_reads):
    '''
    #### Variant type is "INDEL" ####
                    ############################# PLEASE NOTE #######################################
                    ### Variant calling for indels: For insertion, check whether the length of the 
                    ### insertion and sequence matches alt allele. If there is no indel at the anchor
                    ### position (even if the base at the anchor position doesn't match), we consider
                    ### the read as adding support to the reference. For deletions, if the length of
                    ### deletion matches alt allele, we consider that read supporting the alt allele,
                    ### and if there is no deletion at that position (even if there are mismatches in
                    ### the bases spanning the deletion), it's considered to support reference allele.
                    ### Examples for insertions:
                    ### Let's say that the variant is chr1:12345 A > AT
                    ### Scenario1: The read has a C at chr1:12345 along with insertion of T
                    ###            This read will be used to add support to the alternate allele
                    ### Scenario2: Read has a mismatch (let's say 'C') at chr1:12345, but no indel
                    ###            This read will be used to add support to the reference allele
                    ### Scenario3: Read has a different insertion, let's say 'G' insted of 'T'
                    ###            This read will go into the other_reads category.
                    ### Example for deletions:
                    ### Let's say the variant is chr1:12345 AT > A
                    ### Scenario1: The read has a C at chr1:12345 along with deletion of T
                    ###            It will be used to add support to the alt allele.
                    ### Scenario2: The read has a A at chr1:12345 followed by a 2nt deletion
                    ###            This read will go into the other_reads_category
                    ###################################################################################
    '''
    is_alt = False
    if len(ref)==1:
        #Variant is an insertion
        if pileupread.indel == 0: ### and pileupread.alignment.query_sequence[pos_in_read:pos_in_read + 1] == ref:
            ref_reads.append(read_name)
        elif pileupread.indel  == len(alt) - len(ref) and pileupread.alignment.query_sequence[pos_in_read+1:pos_in_read + len(alt)] == alt[1:]:
            alt_reads.append(read_name)
            is_alt = True
            if pileupread.alignment.query_sequence[pos_in_read] != alt[0]:
                anchor_mismatch+=1
        else:
            other_reads.append(read_name)
    else:
        # Variant is a deletion (len(ref)>1 and len(alt)==1)
        if pileupread.indel == 0: ## and pileupread.alignment.query_sequence[pos_in_read+1:pos_in_read + len(ref)] == ref[1:]:
            ref_reads.append(read_name)
        elif pileupread.indel == len(alt) - len(ref): ## and pileupread.alignment.query_sequence[pos_in_read:pos_in_read + len(alt)] == alt:
            alt_reads.append(read_name)
            is_alt = True
            if pileupread.alignment.query_sequence[pos_in_read] != alt[0]:
                anchor_mismatch+=1
        else:
            other_reads.append(read_name)
    return ref_reads, alt_reads, other_reads, anchor_mismatch, is_alt


def compute_vaf(alt_count, dp):
    '''
        Compute VAF from dp and alt_count.
    '''
    if not isinstance(alt_count,int):
        if str(int(alt_count)) == alt_count:
            alt_count=int(alt_count)
        else:
            raise ValueError("alt_count should be an integer. "+alt_count+" provided. Cannot run")
    if not isinstance(dp,int):
        if str(int(dp)) == dp:
            dp=int(dp)
        else:
            raise ValueError("dp should be an integer. "+dp+" provided. Cannot run")
    if alt_count > dp:
        raise ValueError("alt_count {0} is greater than depth {1}.".format(alt_count,dp))
    return (0 if dp==0 else round(float(alt_count)/dp, 4))


def parse_pileup(pileup, pos, ref, alt, variant_type,
                min_mapping_quality=10, min_base_quality=10):
    ''''
        'AD' : 'Allelic depths for the ref and alt alleles in the order listed ', 
        'DP' : 'Depth of coverage: Number of reads covering site computed using a custom pileup-based approach', 
        'AF' : 'Allele fraction of alternate allele in the sample computed using a custom pileup-based approach', 
        'NORMAL_AF' : 'Allele fraction of reference allele in the sample computed using a custom pileup-based approach', 
        'BQ' : 'Base quality for reads supporting the Alt at the start position', 
        'MQ' : 'MAPQ for reads supporting the Alt'
    '''
    bqs = []
    mqs = []
    ref_reads = []
    alt_reads = []
    other_reads = []
    anchor_mismatch=0
    for pileupcolumn in pileup:
        if pileupcolumn.pos == pos - 1:
            for pileupread in pileupcolumn.pileups:
                # if the position in the read is .is_del pos is none so take next
                pos_in_read = pileupread.query_position
                # skip reads where the position is already a deletion (is_del)
                if not pos_in_read:
                    continue
                #  ==========================
                #  pysam filters secondary, dup, and qcfail by default unless nofilter is used
                #  ==========================
                BQ = pileupread.alignment.query_qualities[pos_in_read]
                MQ = pileupread.alignment.mapping_quality
                if BQ < min_base_quality \
                        or MQ < min_mapping_quality \
                        or pileupread.alignment.is_supplementary:
                    continue
                read_name = pileupread.alignment.query_name
                # filter reads that don't span the indel
                if pos_in_read + len(ref) > pileupread.alignment.query_alignment_end \
                        or pos_in_read + len(alt) > pileupread.alignment.query_alignment_end:
                    continue
                is_alt = False
                #  ==========================
                # SNV MNVs
                #  ==========================
                if variant_type == "SNV" or variant_type == "MNV":
                    ref_reads, alt_reads, other_reads, anchor_mismatch, is_alt = parse_snv_mnv(read_name, pileupread, 
                                                                                        ref, alt, pos_in_read, 
                                                                                        anchor_mismatch,
                                                                                        ref_reads, alt_reads, 
                                                                                        other_reads)
                elif variant_type == "COMPLEX":
                    #  ==========================
                    # COMPLEX
                    #  ==========================
                    ### Exact length and sequence of allele match required for complex events. Example: if the variant is AC>T, and if a read has deletion of C but the nt at anchor position is A, it will go into other_reads ###
                    if pileupread.indel == 0 and pileupread.alignment.query_sequence[pos_in_read:pos_in_read + len(ref)] == ref:
                            ref_reads.append(read_name)
                    elif pileupread.indel == len(alt) - len(ref) and pileupread.alignment.query_sequence[pos_in_read:pos_in_read + len(alt)] == alt:
                        alt_reads.append(read_name)
                        is_alt = True
                    else:
                        other_reads.append(read_name)
                else:
                    #  ==========================
                    # INDELs
                    #  ==========================
                    ref_reads, alt_reads, other_reads, anchor_mismatch, is_alt = parse_indel(read_name, pileupread, 
                                                                                        ref, alt, pos_in_read, 
                                                                                        anchor_mismatch,
                                                                                        ref_reads, alt_reads, 
                                                                                        other_reads)
                if is_alt:
                    bqs.append(BQ)
                    mqs.append(MQ)
    ## If there are more than 2 reads that support the alternate allele of an indel variant, but the anchor base does not match, we report that as a PossiblyComplex event ##
    ## Similarly, if there are more than 2 reads that support alt allele of an SNV, but have an indel immediately following the SNV variant, we report that as a PossiblyComplex event ##
    if anchor_mismatch > 2:
        possible_complex=True
    # check sets to make sure reads don't show up in multiple sets
    # supporting multiple calls
    set_ref_raw = set(ref_reads)
    set_alt_raw = set(alt_reads)
    set_other_raw = set(other_reads)
    ref_reads_set = set_ref_raw - set_alt_raw - set_other_raw
    alt_reads_set = set_alt_raw - set_ref_raw - set_other_raw
    other_reads_set = set_other_raw - set_ref_raw - set_alt_raw
    all_reads_set = alt_reads_set|ref_reads_set|other_reads_set
    # tally set in ref and alt, non-ref/alt, all reads
    ref_count = len(ref_reads_set)
    alt_count = len(alt_reads_set)
    # other_count=len(other_reads_set) # not used
    ad = (ref_count, alt_count)
    dp = len(all_reads_set)
    # get VAF
    vaf = compute_vaf(alt_count, dp)
    normal_vaf = compute_vaf(ref_count, dp)
    return ad, dp, vaf, normal_vaf, bqs, mqs

    
    
def annotate(record, normal_bams, non_normal_bams, aliquot_barcodes,
              MAX_INDEL_LEN, min_mapping_quality=10, min_base_quality=10):
    for alt in record.alts:
        too_long = is_too_long(record.ref, alt, MAX_INDEL_LEN)
        if too_long:
            return record
    for aliquot_barcode in aliquot_barcodes:
        if aliquot_barcode in normal_bams:
            bam = normal_bams[aliquot_barcode]
        elif aliquot_barcode in non_normal_bams:
            bam = non_normal_bams[aliquot_barcode]
        pileup = bam.pileup(record.contig, record.pos - 1, 
                             record.pos)
        variant_type=record.info['TYPE']
        ad, dp, vaf, normal_vaf, bqs, mqs = parse_pileup(pileup, 
                                                        record.pos,
                                                        record.ref, 
                                                        record.alts[0],
                                                        variant_type,
                                                        min_mapping_quality, 
                                                        min_base_quality)
        record.samples[aliquot_barcode]['AD'] = ad
        record.samples[aliquot_barcode]['DP'] = dp
        record.samples[aliquot_barcode]['AF'] = vaf
        record.samples[aliquot_barcode]['NORMAL_AF'] = normal_vaf
        record.samples[aliquot_barcode]['AF'] = vaf
        record.samples[aliquot_barcode]['BQ'] = bqs
        record.samples[aliquot_barcode]['MQ'] = mqs
    return record


def count_vcf(bcf_in, bcf_out, 
                normal_bams, non_normal_bams,
                MAX_INDEL_LEN, min_mapping_quality=10, min_base_quality=10):
    '''
        'AD' : 'Allelic depths for the ref and alt alleles in the order listed ', 
        'DP' : 'Depth of coverage: Number of reads covering site computed using a custom pileup-based approach', 
        'AF' : 'Allele fraction of alternate allele in the sample computed using a custom pileup-based approach', 
        'NORMAL_AF' : 'Allele fraction of reference allele in the sample computed using a custom pileup-based approach', 
        'BQ' : 'Base quality for reads supporting the Alt at the start position', 
        'MQ' : 'MAPQ for reads supporting the Alt'
    '''
    aliquot_barcodes =  list(bcf_in.header.samples)
    for record in bcf_in.fetch():
        record = annotate(record, normal_bams, non_normal_bams, 
                          aliquot_barcodes, MAX_INDEL_LEN,
                          min_mapping_quality, min_base_quality)
        exit_status = bcf_out.write(record)
        if exit_status != 0:
            print(exit_status)


def prep_bams(normal_bam_file, non_normal_bam_files):
    normal_bam = read_bam(bam_file=normal_bam_file)
    normal_barcode_aliquot = get_bam_sample(bam=normal_bam)
    normal_bams = {normal_barcode_aliquot : normal_bam}
    non_normal_bams = {}
    for tumor_bam_file in non_normal_bam_files:
        tumor_bam = read_bam(bam_file=tumor_bam_file)
        non_normal_barcode_aliquot = get_bam_sample(bam=tumor_bam)
        non_normal_bams[non_normal_barcode_aliquot] = tumor_bam
    return normal_bams, non_normal_bams


def get_args():
    parser = argparse.ArgumentParser(
                                     description='DESCRIPTION: Runs pileup on tumor(s) \
                                     and normal bam files to compute allele counts for \
                                     bi-allelic SNV and Indel variants in VCF file and \
                                     adds pileup format entries to the VCF file. \
                                     ')
                                     #    Documentation parameters
    parser.add_argument('-v', '--vcf',
                        dest='vcf',
                        help='Input VCF file',
                        required=True)
    parser.add_argument('-t', '--tumor-bams',
                        dest='non_normal_bams',
                        help='Non_normal_bams',
                        nargs='*',
                        required=True)
    parser.add_argument('-n', '--normal-bam',
                        dest='normal_bam',
                        help='Normal bams',
                        required=True)
    parser.add_argument('-o', '--out',
                        dest='out',
                        help='Output VCF file',
                        required=True)
    parser.add_argument('--min-base-quality',
                        dest='min_base_quality',
                        help='Minimum base quality',
                        default=10,
                        required=False)
    parser.add_argument('--min-mapping-quality',
                        dest='min_mapping_quality',
                        help='Minimum mapping quality',
                        default=10,
                        required=False)
    parser.add_argument('--max-indel-len-for-count',
                        dest='max_indel_len_for_count',
                        help='Maximum indel or delin (complex event) length for generating counts',
                        default=10,
                        required=False)
    args = parser.parse_args()
    return args.__dict__
        

def prep_output():
    args = get_args()
    normal_bams, non_normal_bams = prep_bams(normal_bam_file=args['normal_bam'],
                                            non_normal_bam_files=args['non_normal_bams'])
    print(args['vcf'])
    bcf_in = read_vcf(args['vcf'])
    formats = {
            'AD' : 'Allelic depths for the ref and alt alleles in the order listed ', 
            'DP' : 'Depth of coverage: Number of reads covering site computed using a custom pileup-based approach', 
            'AF' : 'Allele fraction of alternate allele in the sample computed using a custom pileup-based approach', 
            'NORMAL_AF' : 'Allele fraction of reference allele in the sample computed using a custom pileup-based approach', 
            'BQ' : 'Base quality for reads supporting the Alt at the start position', 
            'MQ' : 'MAPQ for reads supporting the Alt'
            }
    bcf_in = modify_header(bcf_in)
    bcf_out = pysam.VariantFile(args['out'], 'w', header=bcf_in.header)
    count_vcf(bcf_in, bcf_out, normal_bams, non_normal_bams,
              MAX_INDEL_LEN=args['max_indel_len_for_count'],
              min_mapping_quality=args['min_mapping_quality'], 
              min_base_quality=args['min_base_quality'])
    


def main():
    prep_output()

##########################################################################
#####       Execute main unless script is simply imported     ############
#####                for individual functions                 ############
##########################################################################


if __name__ == '__main__':
    main()
