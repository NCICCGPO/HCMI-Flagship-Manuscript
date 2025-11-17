#!/usr/bin/env python
#    USAGE: python make_pair_vcfs.py PARTICIPANT_FILE OUT_FILE_DIR OUT_FILE_SUFFIX
#   DESCRIPTION:
################################################################################
##################### COPYRIGHT ################################################
# New York Genome Center

# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright (2024) by the New York
# Genome Center. All rights are reserved. This software is supplied without
# any warranty or guaranteed support whatsoever. The New York Genome Center
# cannot be responsible for its use, misuse, or functionality.
# Author: ennifer M Shelton
##################### /COPYRIGHT ###############################################
################################################################################
import sys
import os
import logging as log
import pysam
import argparse
##########################################################################
##############                  Custom functions              ############
##########################################################################


def read_vcf(vcf_file):
    '''
        Read in annotated VCF file.
    '''
    bcf_in = pysam.VariantFile(vcf_file)  # auto-detect input format
    return bcf_in


def test_pass(record):
    '''Test if Mutect2 in Multisample mode filtered out variant'''
    if 'Mutect2Filter' not in record.info.keys():
        return True
    filters = list(record.info['Mutect2Filter'])
    if len(filters) == 1 and filters[0] == 'PASS':
        return True
    return False


def read_support(record, tumors):
    possibles = []
    if 'BQ' in record.format.keys():
        for tumor in tumors:
            if record.samples[tumor]['BQ'][0] != None:
                possible = 'AnyReadSupportsAlt'
            else:
                possible = 'NotAnyReadSupportsAlt'
            possibles.append(possible)
    return all([p == 'AnyReadSupportsAlt' for p in possibles])



def lookup_participant_record(participant_in, record, alt):
    '''
        Check if matching variant is in the GRM VCF.
    '''
    for participant_record in participant_in.fetch(record.contig, record.pos - 1, record.pos):
        if participant_record.ref == record.ref:
            for alt_index, participant_alt in enumerate(participant_record.alts):
                if participant_alt == alt:
                    return participant_record
    return False


def formally_called(record, aliquot_barcode):
    aliquots_called_by = record.info['aliquots_called_by']
    if aliquot_barcode in aliquots_called_by:
        return True
    return False


def keep_var(record, tumors, aliquot_barcode, participant_in):
    '''
        Test if all tumors have > 0 reads (of any BQ)
        which align in all tumor samples.
        Test also that Mutect2 run in Multisample mode 
        forcecalling on the participant VCF passed this variant.
    '''
    participant_record = lookup_participant_record(participant_in, 
                                                    record, 
                                                    record.alts[0])
    mutect2Pass = test_pass(participant_record)
    if not mutect2Pass:
        return False
    formallyCalled = formally_called(record, aliquot_barcode)
    if formallyCalled :
        return True
    readSupportPass = read_support(participant_record, tumors)
    return readSupportPass

def filter_vcf(bcf_in, participant_in, bcf_out, tumors, tumor):
    '''
        Change filter column from PASS to GRM or
        add GRM to filters.
    '''
    for record in bcf_in.fetch():
        keep = keep_var(record, tumors, aliquot_barcode=tumor, participant_in=participant_in)
        if keep:
            exit_status = bcf_out.write(record)
            if exit_status != 0:
                print(exit_status)


def get_args():
    '''Parse input flags
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument('--normal',
                        help='Normal aliquot barcode',
                        required=True
                       )
    parser.add_argument('--participant-vcf',
                        help='Participant VCF',
                        required=True
                       )
    parser.add_argument('--pair-vcf',
                    help='Pair VCF',
                    required=True
                   )
    parser.add_argument('--out-vcf',
                        help='Out VCF file',
                        required=True
               )
    args_namespace = parser.parse_args()
    return args_namespace.__dict__

def main():
    '''
        Change filter column from PASS to GRM or
        add GRM to filters.
    '''
    args = get_args()
    participant_in = read_vcf(args['participant_vcf'])
    bcf_in = read_vcf(args['pair_vcf'])
    samples = list(participant_in.header.samples)
    pair = samples = list(bcf_in.header.samples)
    tumors = [t for t in samples if not t == args['normal']]
    tumor = [t for t in pair if not t == args['normal']][0]
    bcf_out = pysam.VariantFile(args['out_vcf'], 'w', 
                                header=bcf_in.header)
    filter_vcf(bcf_in, participant_in, bcf_out, tumors, tumor)



##########################################################################
#####       Execute main unless script is simply imported     ############
#####                for individual functions                 ############
##########################################################################
if __name__ == '__main__':
    main()
