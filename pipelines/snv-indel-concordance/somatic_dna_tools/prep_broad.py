#!/usr/bin/env python
# USAGE: prep_broad.py VCF_IN VCF_OUT NORMAL TUMOR 
# DESCRIPTION: Print a VCF file with the sample order indicated

import sys
import pandas as pd
import shutil
import logging as log
import os


def load_header(vcf_in):
    '''
        Load a VCF file header as a list of lines.
    '''
    with open(vcf_in) as vcf:
        header = [line for line in vcf if line.startswith('#')]
    return header


def find_header(vcf_file):
        '''
            Get VCF header line numbers (because pandas can't skip
            based on > 1 character words. ## is a comment in VCF but
            # can occur in the VCF INFO fields.
            '''
        try:
            with open(vcf_file) as vcf:
                for i, line in enumerate(vcf):
                    if line.startswith('#'):
                        last = i
        except UnicodeDecodeError:
            with open(vcf_file, encoding='latin-1') as vcf:
                for i, line in enumerate(vcf):
                    if line.startswith('#'):
                        last = i
        return last + 1


def load_vcf(vcf_in, header, normal, tumor):
    '''
        Load a VCF file as an pandas dataframe.
    '''
    last = find_header(vcf_in)
    names = header[-1].rstrip().replace('^#', '').split('\t')
    # 9 10
    if normal.startswith(names[9]):
        names[9] = normal
        names[10] = tumor
    else:
        names[9] = tumor
        names[10] = normal
    vcf_reader = pd.read_csv(vcf_in, skiprows=last,
                             names=names, sep='\t')
    cols = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER',
                 'INFO', 'FORMAT', normal, tumor]
    return vcf_reader[cols]


def vcf_writer(header, vcf_reader, vcf_out_file):
    '''
       Write out the VCF file with corrected sample names
    '''
    with open(vcf_out_file, 'w') as vcf_out:
        for line in header[:-1]:
            vcf_out.write(line)
    vcf_reader.to_csv(vcf_out_file, sep='\t',
                      mode='a', index=False)


def rename(vcf_file, vcf_out_file, normal, tumor):
    '''
        Add prefix to sample name
    '''
    #  =====================
    #  test if renaming should occur
    #  =====================
    rename = False
    if vcf_out_file == vcf_file:
        vcf_out_file = vcf_file + '_tmp.vcf'
        rename = True
    #  =====================
    #  reorder
    #  =====================
    header = load_header(vcf_file)
    vcf_reader = load_vcf(vcf_file, header, normal, tumor)
    #  =====================
    #  rewrite
    #  =====================
    vcf_writer(header, vcf_reader, vcf_out_file)
    #  =====================
    #  rename output VCF
    #  =====================
    if rename:
        shutil.move(vcf_out_file, vcf_file)
    return True


def main():
    vcf_file = sys.argv[1]
    vcf_out_file = sys.argv[2]
    normal = sys.argv[3]
    tumor = sys.argv[4]
    assert os.path.isfile(vcf_file), 'Failed to find prep caller VCF call file :' + vcf_file
    rename(vcf_file, vcf_out_file, normal, tumor)


#  =====================
#  Main
#  =====================


if __name__ == "__main__":
    main()