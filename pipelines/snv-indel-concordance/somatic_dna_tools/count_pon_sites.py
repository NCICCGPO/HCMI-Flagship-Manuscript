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


def count_sites(simple_table_file):
    ''' pon must convert 1-based VCF to 0-based BED
    PON header must include:
    #CHROM, START, END, TotalUniqueSamples
    no other columns are required
    '''
    simple_table = read_table(simple_table_file)
    simple_table['TotalUniqueSamples'] = 1
    simple_table['START'] = simple_table.POS - 1
    simple_table['END'] = simple_table.START + simple_table.ALT.str.len()
    simple_table = simple_table[simple_table.FILTER.isin(['PASS', 'SUPPORT'])].copy()
    simple_table = simple_table[['#CHROM', 'START', 'END', 'TotalUniqueSamples']].drop_duplicates()
    return simple_table


def main():
    simple_table_file = sys.argv[1]
    out_bed = sys.argv[2]
    assert os.path.isfile(simple_table_file), 'Failed to find file :' + simple_table_file
    log.info('Count sites for single')
    pon = count_sites(simple_table_file)
    log.info('Write output file: ' + out_bed)
    pon.to_csv(out_bed, sep='\t', index=None)


if __name__ == "__main__":
    main()