
# USAGE: python calc_baf.py COUNTS_TXT BAF_TXT
# DESCRIPTION: create BAF file for ploting of B-Allele Frequency
import sys
import os
import pandas as pd


def calc_baf(row):
    '''
        Calc BAF from ref_cnt and alt_cnt.
    '''
    if row.T_ref_cnt == 0 and row.T_alt_cnt == 0:
        return 0
    else:
        return row.T_alt_cnt / float(row.T_ref_cnt + row.T_alt_cnt)


def read_counts(file):
    '''
        Read in Allele counts text file
    '''
    data = pd.read_csv(file, sep='\t', dtype={'#Chr' : str})
#    print(data.head())
    data = data[(data.N_dp >= 20)]
    baf_data = data[['#Chr', 'Pos', 'T_ref_cnt', 'T_alt_cnt']].copy()
    if baf_data.empty:
        columns = ['CHR', 'POS', 'REF_COUNT', 'ALT_COUNT', 'BAF']
        baf_data = pd.DataFrame({c:[] for c in columns})
    else:
        baf_data['baf'] = baf_data.apply(lambda row: calc_baf(row), axis=1)
        baf_data.columns = ['CHR', 'POS', 'REF_COUNT', 'ALT_COUNT', 'BAF']
    return baf_data


def write_baf(file, baf_data):
    '''
        Write out a tab delimited file for plots.
    '''
    baf_data.to_csv(file, sep='\t', index=False)


def __main__():
    counts = sys.argv[1]
    baf = sys.argv[2]
    assert os.path.isfile(counts), 'Cannot open file ' + counts
    baf_data = read_counts(counts)
    write_baf(baf, baf_data)

if __name__ == "__main__":
    __main__()
