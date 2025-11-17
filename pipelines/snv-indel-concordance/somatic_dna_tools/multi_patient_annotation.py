#!/usr/bin/env python
#    USAGE: python multi_patient_annotation.py VCF_FILE MULTI_CENTER_OUT_FILE
#   DESCRIPTION: filters out based on how many unique callers their are
################################################################################
##################### COPYRIGHT ################################################
# New York Genome Center

# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright (2024) by the New York
# Genome Center. All rights are reserved. This software is supplied without
# any warranty or guaranteed support whatsoever. The New York Genome Center
# cannot be responsible for its use, misuse, or functionality.
# Version: 0.2
# Author: Jennifer M Shelton
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






def filter_vcf(bcf_in, bcf_filtered_out):
    for record in bcf_in.fetch():
        if record.info['MultiCallerCalled']:
            exit_status = bcf_filtered_out.write(record)
            if exit_status != 0:
                print(exit_status)

def parse_args():
    '''
    '''
    parser = argparse.ArgumentParser(
                                 description='DESCRIPTION: Filter for MultiCallerCalled')
                                 #    Documentation parameters
    parser.add_argument('-v', '--vcf',
                        dest='vcf_file',
                        help='VCF file',
                        required=True)
    parser.add_argument('-f', '--multi-caller-out',
                        dest='multi_caller_out_file',
                        help='Output VCF file for MultiCallerCalled variants',
                        required=True)
    args = parser.parse_args()
    return args

def main():
    '''
       Filter based on MultiCallerCalled.
    '''
    args = parse_args()
    assert os.path.isfile(args.vcf_file), 'Failed to find somatic VCF call file :' + args.vcf_file
    bcf_in = read_vcf(args.vcf_file)
    print(args.multi_caller_out_file)
    bcf_multi_caller_out_file = pysam.VariantFile(args.multi_caller_out_file, 'w',
                                                header=bcf_in.header)
    filter_vcf(bcf_in,bcf_multi_caller_out_file)


##########################################################################
#####       Execute main unless script is simply imported     ############
#####                for individual functions                 ############
##########################################################################
if __name__ == '__main__':
    main()
