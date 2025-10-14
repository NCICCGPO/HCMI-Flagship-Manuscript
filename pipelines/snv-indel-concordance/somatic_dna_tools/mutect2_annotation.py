#!/usr/bin/env python
#	USAGE: python mutect2_annotation.py
#   DESCRIPTION:  Annotates VCF with Mutect2, multi-sample mode forcecalls:
#         1) adds unique FILTER from Mutect2 file to annotated VCF header
#         2) adds INFO from Mutect2 file to annotated VCF header 
#         3) annotates with MUTECT2 GT, FILTER, etc
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


def add_filters_header(mutect2_bcf_in, bcf_out, prefix='Mutect2MultiFilter_'):
    '''
        Add mutect2 filter descriptions (these
        will be in the Mutect2Filter intem in INFO but this keeps the definition
        handy)
    '''
    filters = mutect2_bcf_in.header.filters.keys()
    for filter in filters:
        bcf_out.header.filters.add(id=prefix + filter,
                               number=None,
                               type=None,
                               description=mutect2_bcf_in.header.filters['base_qual'].description)
    return bcf_out


def add_infos_header(mutect2_bcf_in, bcf_out, prefix='Mutect2Multi_'):
    '''
        Add mutect2 INFO descriptions 
    '''
    infos = mutect2_bcf_in.header.info.keys()
    for info in infos:
        number = mutect2_bcf_in.header.info[info].number
        # count of INFO values gets disrupted by multisample mode
        count_disrupted = ['MPOS', 'AS_SB_TABLE', 'AS_FilterStatus',
            'MBQ', 'MFRL', 'MMQ', 'NALOD', 'NLOD', 'POPAF', 'TLOD', 'RPA']
        if info in count_disrupted:
            number='.'
        bcf_out.header.info.add(id=prefix + info,
                    number=number,
                    type=mutect2_bcf_in.header.info[info].type,
                    description=mutect2_bcf_in.header.info[info].description)
    return bcf_out

def add_formats_header(mutect2_bcf_in, bcf_out, prefix='Mutect2Multi_'):
    '''
        Add mutect2 FORMAT descriptions
    '''
    formats = mutect2_bcf_in.header.formats.keys()
    # count of INFO values gets disrupted by multisample mode
    count_disrupted = ['AD', 'AF', 'F1R2', 'F2R1', 'FAD']
    for format in formats:
        number = mutect2_bcf_in.header.formats[format].number
        if format in count_disrupted:
            number='.'
        bcf_out.header.formats.add(id=prefix + format,
                    number=number,
                    type=mutect2_bcf_in.header.formats[format].type,
                    description=mutect2_bcf_in.header.formats[format].description)
    return bcf_out


def match_record(record, mutect2_bcf_in):
    ''' File being annotated must not have multi-allelic rows'''
    assert len(record.alts) == 1, ' File being annotated must not have multi-allelic rows : ' + str(record)
    for mutect2_record in mutect2_bcf_in.fetch(record.contig, 
                                            record.pos - 1, 
                                            record.pos):
        if mutect2_record.ref == record.ref:
            for alt_index, mutect2_alt in enumerate(mutect2_record.alts):
                if mutect2_alt == record.alts[0]:
                    return mutect2_record


def annotate(record, mutect2_bcf_in, prefix='Mutect2Multi_'):
    '''Add Mutect2Multi annotations to record'''
    mutect2_record = match_record(record, mutect2_bcf_in)
    if not mutect2_record:
        # skips MNV splits
        return record
    else:
        # add filters
        record.info['Mutect2Filter'] = mutect2_record.filter.keys()
        # add infos
        for info in mutect2_record.info.keys():
            record.info[prefix + info] = mutect2_record.info[info]
                
        # add formats
        for format in mutect2_record.format.keys():
            for sample in list(mutect2_bcf_in.header.samples):
                if format == 'GT':
                    try:
                        record.samples[sample][format] = mutect2_record.samples[sample][format]
                        record.samples[sample].phased = mutect2_record.samples[sample].phased
                    except ValueError:
                        # Mutect2 called site as multi-allelic 
                        record.info['Mutect2MultiAllelic'] = '|'.join(mutect2_record.alts)
                        record.info['Mutect2MultiAllelicGT'] = '%'.join([str(i) for i in mutect2_record.samples[sample][format]])
                else:
                    record.samples[sample][prefix + format] = mutect2_record.samples[sample][format]
    return record

def annotate_vcf(bcf_in, bcf_out, mutect2_bcf_in):
    '''
    '''
    for record in bcf_in.fetch():
        record = annotate(record, mutect2_bcf_in, prefix='Mutect2Multi_')
        
        exit_status = bcf_out.write(record)
        
        if exit_status != 0:
            print(exit_status)

def main():
    '''
        Annotates VCF with Mutect2, multi-sample mode forcecalls:
        1) adds unique FILTER from Mutect2 file to annotated VCF header
        2) adds INFO from Mutect2 file to annotated VCF header 
        3) annotates with MUTECT2 GT, FILTER, etc
    '''
    #  ==========================
    #  Input variables
    #  ==========================
    parser = argparse.ArgumentParser(
                                     description='DESCRIPTION: Takes in a VCF \
                                     file and preps the file by: \
                                     1) adds unique FILTER from Mutect2 file to annotated VCF header \
                                     2) adds INFO from Mutect2 file to annotated VCF header \
                                     3) annotates with MUTECT2 GT, FILTER, etc \
                                     #    Documentation parameters')
    parser.add_argument('-v', '--vcf',
                        dest='vcf_file',
                        help='VCF file',
                        required=True)
    parser.add_argument('-m', '--mutect2-vcf',
                        dest='mutect2_vcf',
                        help='Tabix indexed Mutect2 vcf file used to annotate vcf',
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
    mutect2_bcf_in = read_vcf(args.mutect2_vcf)
    # Mutect2MultiAllelic
    bcf_in = add_info_header(bcf_out=bcf_in,
                                  id='Mutect2MultiAllelic',
                                  number='1',
                                  type='String',
                                  description='Mutect2 when forcecalling in multi-sample mode called this | sepatated list of alts at this position')
    bcf_in = add_info_header(bcf_out=bcf_in,
                                  id='Mutect2MultiAllelicGT',
                                  number='1',
                                  type='String',
                                  description='Mutect2 GT as string when forcecalling in multi-sample mode called this % sepatated list of alts at this position')

    bcf_in = add_info_header(bcf_out=bcf_in,
                                  id='Mutect2Filter',
                                  number='.',
                                  type='String',
                                  description='Filters applied by Mutect2 when forcecalling in multi-sample mode')
    bcf_in = add_filters_header(mutect2_bcf_in, bcf_in, prefix='Mutect2MultiFilter_')
    bcf_in = add_infos_header(mutect2_bcf_in, bcf_in, prefix='Mutect2Multi_')
    bcf_in = add_formats_header(mutect2_bcf_in, bcf_in, prefix='Mutect2Multi_')
    bcf_out = pysam.VariantFile(args.out, 'w', header=bcf_in.header)
    annotate_vcf(bcf_in, bcf_out, mutect2_bcf_in)


##########################################################################
#####       Execute main unless script is simply imported     ############
#####                for individual functions                 ############
##########################################################################


if __name__ == '__main__':
    main()