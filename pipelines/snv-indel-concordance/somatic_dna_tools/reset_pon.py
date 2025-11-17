#!/usr/bin/env python
#	USAGE: python reset_pon.py
#   DESCRIPTION:
################################################################################
##################### COPYRIGHT ################################################
# New York Genome Center

# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright (2022) by the New York
# Genome Center. All rights are reserved. This software is supplied without
# any warranty or guaranteed support whatsoever. The New York Genome Center
# cannot be responsible for its use, misuse, or functionality.
# Version: 0.1
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


def compose(bcf_in):
    '''
       Remove PON filter.
    '''
    for record in bcf_in.fetch():
        filters = record.filter.keys()
        if 'PON' in filters:
            record.filter.clear()
            if len(filters) == 1:
                record.filter.add('PASS')
            else:
                new_filters = [filter for filter in filters if not filter == 'PON']
                for filter in new_filters:
                    record.filter.add(filter)
        yield record


def write_file(bcf_out, record):
    '''
        Write to a VCF.
    '''
    exit_status = bcf_out.write(record)
    if exit_status != 0:
        print(exit_status)
    return bcf_out


def main():
    '''
        Remove PON filter from a VCF
    '''
    #  ==========================
    #  Input variables
    #  ==========================
    parser = argparse.ArgumentParser(
                                     description='DESCRIPTION: Remove PON filter from a VCF.')
                                     #    Documentation parameters
    parser.add_argument('-v', '--vcf',
                        dest='vcf_file',
                        help='VCF file',
                        required=True)
    parser.add_argument('-o', '--out',
                        dest='out',
                        help='Output VCF file',
                        required=True)
    args = parser.parse_args()
    assert os.path.isfile(args.vcf_file), 'Failed to find caller VCF call file :' + args.vcf_file
    #  ==========================
    #  Undo filter with PON
    #  ==========================
    bcf_in = read_vcf(args.vcf_file)
    bcf_out = pysam.VariantFile(args.out, 'w', header=bcf_in.header)
    bcf_out.header.filters.remove_header('PON')
    for record in compose(bcf_in):
       bcf_out = write_file(bcf_out, record)




##########################################################################
#####       Execute main unless script is simply imported     ############
#####                for individual functions                 ############
##########################################################################


if __name__ == '__main__':
    main()
