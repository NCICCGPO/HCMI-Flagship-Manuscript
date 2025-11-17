#!/usr/bin/env python
#    USAGE: python multi_center_annotation.py VCF_FILE OUT_FILE MULTI_CENTER_OUT_FILE
#   DESCRIPTION: filters out based on how many unique callers their are
################################################################################
##################### COPYRIGHT ################################################
# New York Genome Center

# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright (2023) by the New York
# Genome Center. All rights are reserved. This software is supplied without
# any warranty or guaranteed support whatsoever. The New York Genome Center
# cannot be responsible for its use, misuse, or functionality.
# Version: 0.2
# Author: Kanika Arora (karora@nygenome.org) and Jennifer M Shelton
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


def prep_header(bcf_in, exome_wgs):
    ''' MultiCenterCalled may be adjusted here if NYGC calls are removed but it was added to the header when SNVs were merged back into MNVs (SNVsToMNVs_CountsBasedFilter_AnnotateHighConfCenter.py)
    '''
    bcf_in = add_info_header(bcf_out=bcf_in,
                      id='MultiCallerCalled',
                      number='0',
                      type='Flag',
                      description='Variant called by variant callers (HighConfidence)')
    if exome_wgs:
        bcf_in = add_info_header(bcf_out=bcf_in,
                          id='NotExomeValidated',
                          number='0',
                          type='Flag',
                          description='Variant in exonic region but only called in WGS')
        bcf_in = add_info_header(bcf_out=bcf_in,
                          id='ExomeValidated',
                          number='0',
                          type='Flag',
                          description='Variant in exonic region called in WGS AND Exome (HigherConfidence)')
        bcf_in = add_info_header(bcf_out=bcf_in,
                          id='MissedInWgs',
                          number='0',
                          type='Flag',
                          description='Variant in exonic region but only called by Exome (not the lower coverage WGS alignment)')
        bcf_in = add_info_header(bcf_out=bcf_in,
                          id='IntergenicRegion',
                          number='0',
                          type='Flag',
                          description='Variant is not in exonic region')
    return bcf_in

def annotate(record):
    ''' MultiCenterCalled and MultiCallerCalled'''
    if 'num_centers' in record.info:
        if record.info['num_centers'] > 1:
            record.info['MultiCenterCalled'] = True
        else:
            record.info['MultiCenterCalled'] = False
    if 'called_by' in record.info:
        called_by = list(set(record.info['called_by']))
        if len(called_by) > 1:
            record.info['MultiCallerCalled'] = True
    return record


def annotate_exome_wgs(record):
    ''' NotExomeValidated ExomeValidated
        MissedInWgs IntergenicRegion '''
    exome_called = False
    wgs_called = False
    if record.info['ExonicRegion']:
        if 'Exome_called_by' in record.info:
            exome_count = len(list(set(record.info['Exome_called_by'])))
            if exome_count > 0:
                exome_called = True
        if 'WGS_called_by' in record.info:
            wgs_count = len(list(set(record.info['WGS_called_by'])))
            if wgs_count > 0:
                wgs_called = True
        if exome_called and wgs_called:
            record.info['ExomeValidated'] = True
        elif exome_called:
            record.info['MissedInWgs'] = True
        elif wgs_called:
            record.info['NotExomeValidated'] = True
    else:
        record.info['IntergenicRegion'] = True
    return record


def filter_vcf(bcf_in, bcf_out, bcf_filtered_out, exome_wgs):
    for record in bcf_in.fetch():
        record = annotate(record)
        if exome_wgs:
            record = annotate_exome_wgs(record)
        exit_status = bcf_out.write(record)
        if exit_status != 0:
            print(exit_status)
        if record.info['MultiCallerCalled']:
            exit_status = bcf_filtered_out.write(record)
            if exit_status != 0:
                print(exit_status)

def parse_args():
    '''
        'nygc',
        'nygcWgs',
        'nygcExome',
        'broad',
        'washu',
        'broadWgs',
        'washuWgs',
        'broadExome',
        'washuExome
    '''
    parser = argparse.ArgumentParser(
                                 description='DESCRIPTION: Annotate with relevant details')
                                 #    Documentation parameters
    parser.add_argument('-v', '--vcf',
                        dest='vcf_file',
                        help='VCF file',
                        required=True)
    parser.add_argument('-o', '--out',
                        dest='out',
                        help='Output VCF file',
                        required=True)
    parser.add_argument('-f', '--multi-caller-out',
                        dest='multi_caller_out_file',
                        help='Output VCF file for MultiCallerCalled variants',
                        required=True)
    parser.add_argument('-e', '--exome-wgs',
                        dest='exome_wgs',
                        help='Also annotate if call was called in both Exome and WGS',
                        action='store_true')
    args = parser.parse_args()
    return args

def main():
    '''
        Annotate with MultiCenterCalled and MultiCallerCalled. Filter based on MultiCallerCalled.
    '''
    args = parse_args()
    assert os.path.isfile(args.vcf_file), 'Failed to find somatic VCF call file :' + args.vcf_file
    bcf_in = read_vcf(args.vcf_file)
    bcf_in = prep_header(bcf_in, args.exome_wgs)
    bcf_out = pysam.VariantFile(args.out, 'w', header=bcf_in.header)
    bcf_multi_caller_out_file = pysam.VariantFile(args.multi_caller_out_file, 'w',
                        header=bcf_in.header)
    filter_vcf(bcf_in, bcf_out, bcf_multi_caller_out_file, args.exome_wgs)


##########################################################################
#####       Execute main unless script is simply imported     ############
#####                for individual functions                 ############
##########################################################################
if __name__ == '__main__':
    main()
