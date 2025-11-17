#!/usr/bin/env python
#	USAGE: python consolidate_csq.py --help
#   DESCRIPTION:  Annotates VCF with one CSQ rather than one for each aliquot annotated
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


def determine_csq(bcf_in):
    infos = bcf_in.header.info.keys()
    csqs = [i for i in infos if i.startswith('Aliquot') and i.endswith('_CSQ')]
    descriptions = set([bcf_in.header.info[c].description for c in csqs])
    descriptions = ['Ensembl VEP. Format:' + d.split('Ensembl VEP. Format:')[-1] for d in descriptions]
#     assert len(set(descriptions)) == 1 , 'workflow must merge VCFs with the same annotation'
    if not len(set(descriptions)) == 1:
        if len(csqs) > 0:
#             csqs = [csqs[0]]
            csqs = csqs
        else:
            csqs = []
    if len(csqs) > 0:
        bcf_in.header.info.add(id='CSQ',
                                number=bcf_in.header.info[csqs[0]].number,
                                type=bcf_in.header.info[csqs[0]].type,
                                description=descriptions[0])
    return bcf_in, csqs


def determine_aliquot_map(bcf_in):
    infos = bcf_in.header.info.keys()
    csqs = [i for i in infos if i.startswith('Aliquot') and i.endswith('_CSQ')]
    descriptions = [bcf_in.header.info[c].description for c in csqs]
    aliquot_map = {d.split()[1] : d.split()[0] for d in descriptions}
    print(aliquot_map)
    return aliquot_map

def consolidate(record, csqs):
    for csq in csqs:
        if csq in record.info.keys():
            record.info['CSQ'] = record.info[csq]
    return record

def annotate_vcf(bcf_in, bcf_out, csqs):
    '''
    '''
    for csq in csqs:
        bcf_out.header.info.remove_header(csq)
    for record in bcf_in.fetch():
        record = consolidate(record, csqs)
        for csq in csqs:
            if csq in record.info.keys():
                record.info.__delitem__(csq)
        exit_status = bcf_out.write(record)
        if exit_status != 0:
            print(exit_status)
            
def main():
    '''
        Annotates VCF with one CSQ rather than one for each aliquot annotated
    '''
    #  ==========================
    #  Input variables
    #  ==========================
    parser = argparse.ArgumentParser(
                                     description='DESCRIPTION: Annotates VCF with one CSQ rather than one for each aliquot annotated \
                                     #    Documentation parameters')
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
    bcf_in, csqs = determine_csq(bcf_in)
    aliquot_map = determine_aliquot_map(bcf_in)
    bcf_in.header.add_meta(key='AliquotMap', 
                           value=','.join([':'.join([k, aliquot_map[k]]) for k in aliquot_map]))
    bcf_out = pysam.VariantFile(args.out, 'w', header=bcf_in.header)
    annotate_vcf(bcf_in, bcf_out, csqs)


##########################################################################
#####       Execute main unless script is simply imported     ############
#####                for individual functions                 ############
##########################################################################


if __name__ == '__main__':
    main()
    