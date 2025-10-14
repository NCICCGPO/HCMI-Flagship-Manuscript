#!/usr/bin/env python
#	USAGE: python mutect2_annotation_prep.py
#   DESCRIPTION: Resets FILTER to PASS 
#       and adds FILTER value to the INFO
#       also renames VCF INFO and FORMAT to keep unique
################################################################################
##################### COPYRIGHT ################################################
# New York Genome Center

# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright (2024) by the New York
# Genome Center. All rights are reserved. This software is supplied without
# any warranty or guaranteed support whatsoever. The New York Genome Center
# cannot be responsible for its use, misuse, or functionality.
# Version: 0.1
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


def add_filter_header(bcf_out,
                      id,
                      description):
    '''
        Add new FILTER field
    '''
    bcf_out.header.filters.add(id=id,
                               number=None,
                               type=None,
                               description=description)
    return bcf_out


def add_filters_header(bcf_out, prefix=''):
    '''
        Remove all filters except PASS
    '''
    for id in bcf_out.header.filters.keys():
        if not id in ['PASS']:
            bcf_out.header.filters.remove_header(id)
    return bcf_out


def main():
    '''
        Prepare the VCF to be used in annotation:
        1) adds FILTER value to the INFO
        2) Resets FILTER to PASS 
        3) renames VCF INFO and FORMAT to keep unique
    '''
    #  ==========================
    #  Input variables
    #  ==========================
    parser = argparse.ArgumentParser(
                                     description='DESCRIPTION: Takes in a VCF \
                                     file and preps the file by: \
                                     1) adds FILTER value to the INFO \
                                     2) Resets FILTER to PASS  \
                                     3) renames VCF INFO and FORMAT to keep unique \
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
    #  Run prep
    #  ==========================
    bcf_in = read_vcf(args.vcf_file)
    bcf_in = add_info_header(bcf_out=bcf_in,
                                  id='mutect2Filter',
                                  number='.',
                                  type='String',
                                  description='Filters applied by Mutect2 when forcecalling in multi-sample mode')
                                  

    bcf_out = pysam.VariantFile(args.out, 'w', header=bcf_in.header)
    if filter:
        bcf_out = remove_filters_header(bcf_out)
    write_file(bcf_in,
               bcf_out,
               tool=args.tool,
               filter=filter,
               skip_called_by=args.skip_called_by,
               exome_wgs=args.exome_wgs)


##########################################################################
#####       Execute main unless script is simply imported     ############
#####                for individual functions                 ############
##########################################################################


if __name__ == '__main__':
    main()