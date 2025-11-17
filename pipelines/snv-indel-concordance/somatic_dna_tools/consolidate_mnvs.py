#!/usr/bin/env python
#	USAGE: python consolidate_mnvs.py
#   DESCRIPTION: Takes in a VCF \
#          file and preps the file by: \
#          1) setting v6 strelka2 annotation to match centers \
#          2) Deduplicate SNVs from MNVs
#          3) Reset MNVs and split MNVs to PASS
# ##############################################################################
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
called_by_pattern = re.compile('Aliquot\d+_called_by')

def read_vcf(vcf_file):
    '''
        Read in annotated VCF file.
        '''
    bcf_in = pysam.VariantFile(vcf_file)  # auto-detect input format
    return bcf_in

def get_mergable(bcf_in):
    'send mnv_id if id has same called_by for all records'
    mergable_ids = {}
    infos = bcf_in.header.info.keys()
    aliquots = [i for i in infos if called_by_pattern.match(i)]
    for record in bcf_in.fetch():
        if 'MNV_ID' in record.info.keys():
            for mnv_id in record.info['MNV_ID']:
                if mnv_id not in mergable_ids:
                    mergable_ids[mnv_id] = { 'index': [],
                                             'called_by' : {},
                                             'type' : [],
                                             'variant' : str(record)
                                        }
                    mergable_ids[mnv_id]['index'].append(0)
                else:
                    mergable_ids[mnv_id]['index'].append(mergable_ids[mnv_id]['index'][-1] + 1)
                for aliquot in aliquots:
                    if not aliquot in mergable_ids[mnv_id]['called_by']:
                            mergable_ids[mnv_id]['called_by'][aliquot] = []
                    if aliquot in record.info.keys():
                        mergable_ids[mnv_id]['called_by'][aliquot].append(record.info[aliquot])
                    else:
                        mergable_ids[mnv_id]['called_by'][aliquot].append(tuple([]))
                mergable_ids[mnv_id]['type'].append(record.info['TYPE'])
    return mergable_ids

def filter_mergable(mergable_ids):
    ''' 
        Confirm that MNV has exact same callers as each SNV (othewise print all records)
    '''
    not_mergables = []
    mergables = []
    for mnv_id in mergable_ids:
        print(mnv_id)
        print(type(mnv_id), type(mnv_id) == tuple)
        print(mergable_ids[mnv_id])
        mnv_index = [index for index, type in enumerate(mergable_ids[mnv_id]['type']) if type == 'MNV'][0]
        for aliquot in mergable_ids[mnv_id]['called_by']:
            mnv_called_by = mergable_ids[mnv_id]['called_by'][aliquot][mnv_index]
            snv_called_bys = [called_by for index, called_by in enumerate(mergable_ids[mnv_id]['called_by'][aliquot]) if index != mnv_index][0]
            # test if any snv has a different number of callers from the mnv for this aliquot
            if len(set(mnv_called_by).symmetric_difference(set(snv_called_bys))) > 0:
                not_mergables.append(mnv_id)
        if not mnv_id in not_mergables:
            mergables.append(mnv_id)
    return not_mergables, mergables
        

def consolidate(bcf_in, bcf_out):
    mergable_ids = get_mergable(bcf_in)
    not_mergables, mergables = filter_mergable(mergable_ids)
    for record in bcf_in.fetch():
        write = True
        if 'MNV_ID' in record.info.keys():
            for mnv_id in record.info['MNV_ID']:
                if (mnv_id in mergables) & (record.info['TYPE'] == 'SNV'):
                    write = False
        if write:
            MultiCallerCalled = False
            record.filter.clear()
            record.filter.add('PASS')
            aliquots = [i for i in record.info.keys() if called_by_pattern.match(i)]
            for aliquot in aliquots:
                record.info[aliquot] = [caller.replace('strelka2_SNVs', 'strelka2') for caller in record.info[aliquot]]
                if len(set(record.info[aliquot])) > 1:
                    MultiCallerCalled = True
            record.info['MultiCallerCalled'] = MultiCallerCalled
            exit_status = bcf_out.write(record)
            if exit_status != 0:
                print(exit_status)

def main():
    '''
        DESCRIPTION: Takes in a VCF \
         file and preps the file by: \
         1) setting v6 strelka2 annotation to match centers \
         2) Deduplicate SNVs from MNVs
         3) Reset MNVs and split MNVs to PASS
    '''
    #  ==========================
    #  Input variables
    #  ==========================
    parser = argparse.ArgumentParser(
                                     description='DESCRIPTION: Takes in a VCF \
                                     file and preps the file by: \
                                     1) setting v6 strelka2 annotation to match centers \
                                     2) Deduplicate SNVs from MNVs \
                                     3) Reset MNVs and split MNVs to PASS \
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
    bcf_in.header.info.add(id='MultiCallerCalled',
                            number='0',
                            type='Flag',
                            description='Variant called by variant callers in at least one aliquot (HighConfidence)')
    bcf_out = pysam.VariantFile(args.out, 'w', header=bcf_in.header)
    consolidate(bcf_in, bcf_out)


##########################################################################
#####       Execute main unless script is simply imported     ############
#####                for individual functions                 ############
##########################################################################


if __name__ == '__main__':
    main()
    