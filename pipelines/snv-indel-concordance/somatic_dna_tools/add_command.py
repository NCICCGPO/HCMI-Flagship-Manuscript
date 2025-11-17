#!/usr/bin/env python

# env:
#  - python >=3.6
#  - pysam >= 0.11.2.1
# run_with: python

# USAGE: sample_check.py VCF [TUMOR] NORMAL
# DESCRIPTION: adds commands and a timestamp to VCF files. Must use pysam
# 0.11.2.1 if using PINDEL because later versions consider "END" reserved and therefore
# will ommit this when printing results from Pindel.

import sys
import time
import datetime
import json
import shutil
import pysam
from dateutil.parser import parse


def add_command(vcf_file, vcf_out_file, line_dict):
    '''
        Add comment lines (with ID) to a VCF file.
        Line must not include | as these are used
        to split the line into a dictionary.
    '''
    #  =====================
    #  test if renaming should occur
    #  =====================
    rename = False
    if vcf_out_file == vcf_file:
        vcf_out_file = vcf_file + '_tmp.vcf'
        rename = True
    #  =====================
    #  Add command lines
    #  =====================
    bcf_in = pysam.VariantFile(vcf_file)  # auto-detect input format
    cmdlines = []
    for key in line_dict:
        try:
            if parse(line_dict[key]):
                datetime = line_dict[key]
            else:
                cmdlines.append(line_dict[key].replace('"', '\''))
        except ValueError:
            cmdlines.append(line_dict[key].replace('"', '\''))
    for cmdline in cmdlines:
        bcf_in.header.add_meta(key='cmdline', value='<CMD="' + cmdline + '",DATETIME="' + datetime + '">')
    #  =====================
    #  Write out lines
    #  =====================
    bcf_out = pysam.VariantFile(vcf_out_file, 'w', header=bcf_in.header)
    for record in bcf_in.fetch():
        exit_status = bcf_out.write(record)
        if exit_status != 0:
            print(exit_status)
    #  =====================
    #  rename output VCF
    #  =====================
    if rename:
        vcf_test_file = vcf_file + '_pre_command.vcf'
        shutil.copy(vcf_file, vcf_test_file) # temp copy to compare
        shutil.move(vcf_out_file, vcf_file)


def read_lines(json_file):
    '''
        Read JSON format command lines.
    '''
    lines = json.load(open(json_file))
    return lines


def compose_lines(lines):
    '''
        Add lines to indicate when job ran.
    '''
    time_stamp = datetime.datetime.now().strftime("%d/%m/%y %H:%M")
    lines['calling_variants_date'] = time_stamp
    return lines


def add_lines(vcf_file, vcf_out_file, json_file):
    '''
        Add info on commands run to VCF for strelka.
    '''
    lines = read_lines(json_file)
    lines = compose_lines(lines)
    add_command(vcf_file, vcf_out_file, lines)
    return True


#  =====================
#  Main
#  =====================
if __name__ == "__main__":
    vcf_file = sys.argv[1]
    vcf_out_file = sys.argv[2]
    json_file = sys.argv[3]
    add_lines(vcf_file, vcf_out_file, json_file)
