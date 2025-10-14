#!/usr/bin/env python
#	USAGE: python cancer_gene_census.py cancer_gene_census.csv VCF VCF_OUT
#   DESCRIPTION: Annotates files by adding information about the
# Cosmic Genome Census entry for the nearest gene.
################################################################################
##################### COPYRIGHT ################################################
# New York Genome Center

# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright (2018) by the New York
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
import pprint
from collections import OrderedDict
##########################################################################
##############                  Custom functions              ############
##########################################################################


def remove_info(bcf_in):
    '''
        Remove a INFO field from VCF.
    '''
    for id in bcf_in.header.info.keys():
        if id in ['washu_CSQ', 'nygc_CSQ', 'broad_CSQ']:
            bcf_in.header.info.remove_header(id)
    return bcf_in


def remove_format(bcf_in):
    '''
        Remove a FORMAT field from VCF.
    '''
    for id in bcf_in.header.formats.keys():
        if id in ['nygc_AD','nygc_DP', 'nygc_AF']:
            bcf_in.header.formats.remove_header(id)
    return bcf_in


class Variant(object):
    
    
    def __init__(self, record):
        self.record = record
        self.line = str(self.record).rstrip()
        self.parts = self.line.split('\t')
        # VCF columns
        self.chrom = self.parts[0]
        self.pos = self.parts[1]
        self.id = self.parts[2]
        self.ref = self.parts[3]
        self.alts = self.parts[4].split(',')
        self.qual = self.parts[5]
        self.filters = self.parts[6].split(';')
        self.info = self.parts[7].split(';')
        self.format = self.parts[8].split(':')
        self.samples = self.parts[9:]
        # modify
        self.bad_format = ['nygc_AD','nygc_DP', 'nygc_AF']
        self.bad_info = ['washu_CSQ', 'nygc_CSQ', 'broad_CSQ']
        self.info_dict = self.get_info()
        self.samples[0] = self.fix_format(self.samples[0].split(':'))
        self.samples[1] = self.fix_format(self.samples[1].split(':'))
        self.format = [key for key in self.format if not key in self.bad_format]
    
    def get_info(self):
        '''
            Get current info line and add prefix if needed
            '''
        info_dict = OrderedDict()
        for item in self.info:
            if item.split('=')[0] in self.good_info:
                if '=' in item:
                    info_dict.update({item.split('=')[0] : item.split('=')[1]})
                else:
                    info_dict.update({item : None })
        return info_dict
    
    def fix_format(self, sample):
        '''
            Reduce to good formats
        '''
        format_dict = dict(zip(self.format, sample))
        new_format = [format_dict[key] for key in format_dict if not key in self.bad_format]
        return ':'.join(new_format)

    def write(self):
        if ':'.join(self.format) == '':
            self.format = '.'
            self.samples = ['.', '.']
        if ';'.join(self.filters) == 'PASS':
            line = [self.chrom,
                    self.pos,
                    self.id,
                    self.ref,
                    ','.join(self.alts),
                    str(self.qual),
                    ';'.join(self.filters),
                    ';'.join(['='.join([x for x in [key, self.info_dict[key]] if x != None]) for key in self.info_dict]),
                    ':'.join(self.format)]
            line += self.samples
            self.new_line = '\t'.join(line)
            return self.new_line
        else:
            return False


def read_vcf(vcf_file):
    '''
        Read in annotated VCF file.
    '''
    bcf_in = pysam.VariantFile(vcf_file)  # auto-detect input format
    return bcf_in


def write_vcf(bcf_in, vcf_out_file:
    '''
        Write out the download
    '''
    # Import the header after removal of extra metadata
    header = str(bcf_in.header).rstrip()
    # Write new header with fewer metadata keys overall
    with open(vcf_out_file, 'w') as vcf_out:
        for line in header.split('\n'):
            vcf_out.write(line + '\n')
        for record in bcf_in:
            line = Variant(record).write()
            if line:
                vcf_out.write(line + '\n')


def main():
    '''
        Reduce metadata in VCF for main VCF output
    '''
    vcf_file = sys.argv[1]
    vcf_out_file = sys.argv[2]
    bcf_in = read_vcf(vcf_file)
    bcf_in = remove_format(bcf_in)
    bcf_in = remove_info(bcf_in)
    write_vcf(bcf_in, vcf_out_file)


##########################################################################
#####       Execute main unless script is simply imported     ############
#####                for individual functions                 ############
##########################################################################
if __name__ == '__main__':
    main()
