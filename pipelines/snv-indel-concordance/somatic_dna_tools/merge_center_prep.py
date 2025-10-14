#!/usr/bin/env python
#	USAGE: python merge_prep.py
#   DESCRIPTION: Renames VCF INFO and FORMAT to keep unique
################################################################################
##################### COPYRIGHT ################################################
# New York Genome Center

# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright (2018) by the New York
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


def remove_filters_header(bcf_out):
    '''
        Remove all filters except PASS
    '''
    for id in bcf_out.header.filters.keys():
        if not id in ['PASS']:
            bcf_out.header.filters.remove_header(id)
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

def get_nygc_callers(record, exome_wgs=False):
    exome_callers = []
    wgs_callers = []
    callers = []
    if exome_wgs:
        if 'nygc_Exome_called_by' in record.info:
            exome_callers = record.info['nygc_Exome_called_by']
        if 'nygc_WGS_called_by' in record.info:
            wgs_callers = record.info['nygc_WGS_called_by']
            callers += wgs_callers
        if 'nygc_WGS_supported_by' in record.info:
            wgs_callers += record.info['nygc_WGS_supported_by']
        if 'nygc_Exome_supported_by' in record.info:
            exome_callers += record.info['nygc_Exome_supported_by']
        if exome_callers:
            callers += exome_callers
        if wgs_callers:
            callers += wgs_callers
    else:
        callers = record.info['nygc_called_by']
        if 'nygc_supported_by' in record.info:
            callers += record.info['nygc_supported_by']
    results = {
        'callers' : callers,
        'exome_callers' : exome_callers,
        'wgs_callers' : wgs_callers
        }
    return results

def get_broad_callers(record, exome_wgs=False):
    exome_callers = False
    wgs_callers = False
    callers = []
    print(record)
    assert len(record.alts) == 1, 'Vcf must contain one variant per line for this annotation'
    if len(record.alts[0]) == 1 and len(record.ref) == 1:
        callers = ('mutect',)
    else:
        callers = ('strelka2',)
    return callers


def prep_record(record, tool, passing, skip_called_by, exome_wgs=False):
    '''
        Pass lines that have no special characters in ALT/REF
    '''
    record.id = None
    record.qual = None
    if passing and not skip_called_by:
        if 'nygc' in tool:
            nygc_results = get_nygc_callers(record, exome_wgs=exome_wgs)
            callers = list(set(nygc_results['callers']))
        elif 'broad' in tool:
            callers = get_broad_callers(record, exome_wgs=exome_wgs)
        elif 'washu' in tool:
            if 'Exome' in tool:
                callers = record.info['washuExome_set'].split('-')
            elif 'Wgs' in tool:
                callers = record.info['washuWgs_set'].split('-')
            else:
                callers = record.info['washu_set'].split('-')
            callers = [c.replace('strelka', 'strelka2') for c in callers]
        record.info['centers_called_by'] = (tool,)
        record.info['called_by'] = callers
        record.info['num_centers'] = 1
        if exome_wgs:
            if 'washu' or 'broad' in tool:
                if 'Exome' in tool:
                    record.info['Exome_called_by'] = callers
                else:
                    record.info['WGS_called_by'] = callers
            else:
                record.info['Exome_called_by'] = nygc_results['exome_callers']
                record.info['WGS_called_by'] = nygc_results['wgs_callers']
    return record


def write_file(bcf_in, bcf_out, tool, filter=True, skip_called_by=False, exome_wgs=False):
    '''
       Filter based on FILTER column,
       also filter lines with special characters in ALT/REF
    '''
    passing = False
    for record in bcf_in.fetch():
        filters = record.filter.keys()
        #  ====================
        #  Passing variants
        #  ====================
        if (len(filters) == 1 and \
                filters[0] == 'PASS') \
                or len(filters) == 0:
            passing = True
            write = True
        #  ====================
        #  Failing variants
        #  ====================
        else:
            if filter:
                write = False
            else:
                write = True
        if write and not pass_alleles(record):
            write = False
        if write:
            record = prep_record(record, tool, passing, skip_called_by, exome_wgs)
            exit_status = bcf_out.write(record)
            if exit_status != 0:
                print(exit_status)
    return True

def add_headers(bcf_in, prefix=''):
        bcf_in = add_info_header(bcf_out=bcf_in,
                                  id=prefix + 'called_by',
                                  number='.',
                                  type='String',
                                  description='Name of the variant caller(s) that the variant was called by')
        return bcf_in

def main():
    '''
        Prepare the VCF file for merging by:
        1) 'TYPE' is added to header
        2) 'centers_called_by' is added to header
        3) 'num_centers' is added to header
        4) filter lines with special characters in REF/ALT (e.g. <DEL> in manta)
        5) fill in centers_called_by, num_centers, center + '_called_by'
        7) non 'PASS' FILTER lines are removed (if not skip-filter)
    '''
    #  ==========================
    #  Input variables
    #  ==========================
    parser = argparse.ArgumentParser(
                                     description='DESCRIPTION: Takes in a VCF \
                                     file and preps the file by: \
                                     1) "TYPE" is added to header \
                                     2) "centers_called_by" is added to header \
                                     3) "num_centers" is added to header \
                                     4) filter lines with special characters in REF/ALT (e.g. <DEL> in manta) \
                                     5) fill in centers_called_by, num_centers, center + _called_by \
                                     7) non "PASS" FILTER lines are removed (if not skip-filter) \
                                     . Command-line \
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
    parser.add_argument('-t', '--tool',
                        dest='tool',
                        choices=['nygc',
                                 'nygcWgs',
                                 'nygcExome',
                                 'broad',
                                 'washu',
                                 'broadWgs',
                                 'washuWgs',
                                 'broadExome',
                                 'washuExome'],
                        help='Tool name',
                        required=True)
    parser.add_argument('-e', '--exome-wgs',
                        dest='exome_wgs',
                        help='Also annotate if call was called in both Exome and WGS',
                        action='store_true')
    parser.add_argument('-f', '--skip-filter',
                        dest='skip_filter',
                        help='Remove calls that are not PASS',
                        action='store_true')
    parser.add_argument('-s', '--skip-called-by',
                        dest='skip_called_by',
                        help='Skip noting the current center name or call count.',
                        action='store_true')
    args = parser.parse_args()
    filter = True
    if args.skip_filter:
        filter = False
    assert os.path.isfile(args.vcf_file), 'Failed to find caller VCF call file :' + args.vcf_file
    #  ==========================
    #  Run prep
    #  ==========================
    bcf_in = read_vcf(args.vcf_file)
    if not args.skip_called_by:
        bcf_in = add_headers(bcf_in, prefix='')
        if args.exome_wgs:
            bcf_in = add_headers(bcf_in, prefix='Exome_')
            bcf_in = add_headers(bcf_in, prefix='WGS_')
        bcf_in = add_info_header(bcf_out=bcf_in,
                          id='centers_called_by',
                          number='.',
                          type='String',
                          description='Name of the center(s) that the variant was called by')
        bcf_in = add_info_header(bcf_out=bcf_in,
                                  id='num_centers',
                                  number='1',
                                  type='Integer',
                                  description='Number of centers that called the variant')
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
