#!/usr/bin/env python
#    USAGE: python merge_pon_sites.py
#   DESCRIPTION: merge single sample PON files
################################################################################
##################### COPYRIGHT ################################################
# New York Genome Center

# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright (2022) by the New York
# Genome Center. All rights are reserved. This software is supplied without
# any warranty or guaranteed support whatsoever. The New York Genome Center
# cannot be responsible for its use, misuse, or functionality.
# Version: 0.1
# Author: Jennifer M Shelton
##################### /COPYRIGHT ###############################################
################################################################################
import pandas as pd
import sys
import logging as log
import os


def read_table(simple_table_file):
    simple_table = pd.read_csv(simple_table_file, sep='\t')
    return simple_table


def merge_pons(single_pons):
    ''' pon must convert 1-based VCF to 0-based BED
    PON header must include:
    #CHROM, START, END, TotalUniqueSamples
    no other columns are required
    '''
    pon = read_table(single_pons[0])
    for single_pon_file in single_pons[1:]:
        single_pon = read_table(single_pon_file)
        pon = pd.merge(pon, single_pon, on=['#CHROM', 'START', 'END'], how='outer')
        pon = pon.fillna({'TotalUniqueSamples_x' : 0, 'TotalUniqueSamples_y' : 0})
        pon['TotalUniqueSamples'] = pon['TotalUniqueSamples_x'] + pon['TotalUniqueSamples_y']
        pon = pon[['#CHROM', 'START', 'END', 'TotalUniqueSamples']].copy()
    return pon

def main():
    out_bed = sys.argv[1]
    single_pons = sys.argv[2:]
    for pon_file in single_pons:
        assert os.path.isfile(pon_file), 'Failed to find input PON file :' + pon_file
    pon = merge_pons(single_pons)
    pon.to_csv(out_bed, sep='\t', index=None)


if __name__ == "__main__":
    main()