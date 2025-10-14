#!/usr/bin/env python
#	USAGE: python merge_aliquot_prep.py --help
#   DESCRIPTION: Prepare the VCF file for merging by:
#          1) "aliquots_called_by" is added to header \
#          2) "non_normal_kinds_called_by" is added to header \
#          3) fill in aliquots_called_by, non_normal_kinds_called_by \
################################################################################
##################### COPYRIGHT ################################################
# New York Genome Center

# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright (2023) by the New York
# Genome Center. All rights are reserved. This software is supplied without
# any warranty or guaranteed support whatsoever. The New York Genome Center
# cannot be responsible for its use, misuse, or functionality.
# Version: 0.1
# Author: Kanika Arora (karora@nygenome.org) and Jennifer M Shelton
##################### /COPYRIGHT ###############################################
################################################################################
import sys
import os
import logging as log
import pysam
import pandas as pd
import re
import argparse
##########################################################################
##############                  Custom functions              ############
##########################################################################
base_pattern = re.compile(r'^[ACGTN]+$')



def read_vcf(vcf_file):
    '''
        Read in annotated VCF file.
        '''
    bcf_in = pysam.VariantFile(vcf_file)  # auto-detect input format
    return bcf_in


def add_info_header(bcf_out,
                    id,
                    number,
                    type,
                    description):
    '''
        Add new INFO field
    '''
    bcf_out.header.info.add(id=id,
                            number=number,
                            type=type,
                            description=description)
    return bcf_out


def pass_alleles(record, base_pattern=base_pattern):
    '''
        Pass lines that have no special characters in REF/ALT
    '''
    passed = True
    alleles = list(record.alts) + [record.ref]
    for allele in alleles:
        if not re.match(base_pattern, allele):
            passed = False
    return passed


def prep_record(record, aliquot, non_normal_kind):
    '''
        Pass lines that have no special characters in ALT/REF
    '''
    record.info['aliquots_called_by'] = (aliquot,)
    record.info['non_normal_kinds_called_by'] = (non_normal_kind,)
    return record


def write_file(bcf_in, bcf_out, aliquot, non_normal_kind):
    '''
    Write lines with custom annotation
    '''
    passing = False
    for record in bcf_in.fetch():
        record = prep_record(record, aliquot, non_normal_kind)
        exit_status = bcf_out.write(record)
        if exit_status != 0:
            print(exit_status)
    return True


def main():
    '''
        Prepare the VCF file for merging by:
         1) "aliquots_called_by" is added to header \
         2) "non_normal_kinds_called_by" is added to header \
         3) fill in aliquots_called_by, non_normal_kinds_called_by \
    '''
    #  ==========================
    #  Input variables
    #  ==========================
    parser = argparse.ArgumentParser(
                                     description='DESCRIPTION: Takes in a VCF \
                                     file and preps the file by: \
                                     1) "aliquots_called_by" is added to header \
                                     2) "non_normal_kinds_called_by" is added to header \
                                     3) fill in aliquots_called_by, non_normal_kinds_called_by \
                                     options that may be omitted (i.e. are NOT \
                                     required) are shown in square brackets.')
                                     #    Documentation parameters
    parser.add_argument('-v', '--vcf',
                        dest='vcf_file',
                        help='VCF file',
                        required=True)
    parser.add_argument('-o', '--out',
                        dest='out',
                        help='Output VCF file',
                        required=True)
    parser.add_argument('-t', '--aliquot',
                        dest='aliquot',
                        help='Aliquot name',
                        required=True)
    parser.add_argument('-n', '--non-normal-kind',
                        dest='non_normal_kind',
                        help='Non-normal kind',
                        required=True)
    args = parser.parse_args()
    assert os.path.isfile(args.vcf_file), 'Failed to find caller VCF call file :' + args.vcf_file
    #  ==========================
    #  Run prep
    #  ==========================
    bcf_in = read_vcf(args.vcf_file)
    bcf_in = add_info_header(bcf_out=bcf_in,
                      id='aliquots_called_by',
                      number='.',
                      type='String',
                      description='Name of the non-normal aliquot(s) that the variant was called by')
    bcf_in = add_info_header(bcf_out=bcf_in,
                              id='non_normal_kinds_called_by',
                              number='.',
                              type='String',
                              description='The type of non-normal aliquot(s) that the variant was called by (e.g. Tumor, Model)')
    bcf_out = pysam.VariantFile(args.out, 'w', header=bcf_in.header)
    write_file(bcf_in,
               bcf_out,
               aliquot=args.aliquot,
               non_normal_kind=args.non_normal_kind)


##########################################################################
#####       Execute main unless script is simply imported     ############
#####                for individual functions                 ############
##########################################################################


if __name__ == '__main__':
    main()
