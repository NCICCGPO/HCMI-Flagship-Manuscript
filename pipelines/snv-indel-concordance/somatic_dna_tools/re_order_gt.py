#!/usr/bin/env python
#	USAGE: python merge_prep.py
#   DESCRIPTION:
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



##########################################################################
##############                  Custom functions              ############
##########################################################################



class Variant():
    '''
        Import a pysam record. and write out from record.
        The class allows editing of fixed elements like
        the min number of samples in the VCF.
    '''

    def __init__(self, line):
        self.line = line.strip()
        self.parts = self.line.split('\t')
        # VCF columns
        self.chrom = self.parts[0]
        self.pos = self.parts[1]
        self.id = self.parts[2]
        self.ref = self.parts[3]
        self.alts = self.parts[4]
        self.qual = self.parts[5]
        self.filters = self.parts[6]
        self.info = self.parts[7]
        self.format_keys =  self.parts[8].split(':')
        self.samples = self.parts[9:]
        self.set_orders()

    def set_sample_order(self, sample, 
                        gt_like, ordered_keys):
        values = sample.split(':')
        map = {self.format_keys[i]: values[i] for i in range(len(values))}
        for gt in gt_like:
            if map[gt] == '.':
                map[gt] = './.'
        return ':'.join([map[k] for k in ordered_keys])
    
    def set_orders(self):
        ''' set key order (GT and GT like first)'''
        gt_like = [key for key in self.format_keys if key.endswith('GT')]
        ordered_keys = self.order_keys(self.format_keys)
        self.format = ':'.join(ordered_keys)
        if not self.format:
            self.format = '.'
        self.reordered_samples = []
        for sample in  self.samples:
            reordered_sample = self.set_sample_order(sample, gt_like, ordered_keys)
            self.reordered_samples.append(reordered_sample)

    def order_keys(self, format_keys):
        '''
            Return list of keys with values for FORMAT
        '''
        ordered_keys = [key for key in format_keys if key == 'GT']
        ordered_keys += [key for key in format_keys if key.endswith('GT') and key != 'GT']
        ordered_keys += [key for key in format_keys if not key.endswith('GT')]
        return ordered_keys

    def write(self):
        '''
            Return a reformatted string from a pysam object.
        '''
        line = [self.chrom,
                self.pos,
                self.id,
                self.ref,
                self.alts,
                self.qual,
                self.filters,
                self.info,
                self.format]
        line += self.reordered_samples
        self.new_line = '\t'.join(line) + '\n'
        return self.new_line



def load_vcf(vcf_in):
    '''
        Load a VCF file and fixes lines.
    '''
    with open(vcf_in) as vcf:
        for line in vcf:
            if not line.startswith('#'):
                yield Variant(line).write()
            else:
                yield line


def vcf_writer(vcf_reader, vcf_out_file):
    '''
       Write out a VCF file with reordered FORMAT keys
    '''
    with open(vcf_out_file, 'w') as vcf_out:
        for line in vcf_reader:
            vcf_out.write(line)
    return True


def re_order(vcf_file, vcf_out_file):
    '''
        Reorder FORMAT keys
    '''
    #  =====================
    #  test if renaming should occur
    #  =====================
    rename = False
    if vcf_out_file == vcf_file:
        vcf_out_file = vcf_file + '_tmp.vcf'
        rename = True
    #  =====================
    #  rename
    #  =====================
    vcf_reader = load_vcf(vcf_file)
    #  =====================
    #  rewrite
    #  =====================
    vcf_writer(vcf_reader, vcf_out_file)
    #  =====================
    #  rename output VCF
    #  =====================
    if rename:
        shutil.move(vcf_out_file, vcf_file)
    return True



def main():
    '''
    '''
    #  ==========================
    #  Input variables
    #  ==========================
    vcf_file= sys.argv[1]
    vcf_out_file = sys.argv[2]

    assert os.path.isfile(vcf_file), 'Failed to find caller VCF call file :' + vcf_file
    #  ==========================
    #  Run prep
    #  ==========================
    re_order(vcf_file, vcf_out_file)


##########################################################################
#####       Execute main unless script is simply imported     ############
#####                for individual functions                 ############
##########################################################################


if __name__ == '__main__':
    main()