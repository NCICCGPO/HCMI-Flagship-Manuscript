#!/usr/bin/env python
#	USAGE: python reset_mantis.py MANTIS OUT NORMAL
#   DESCRIPTION: Reset Mantis status based on threshhold for unmatched tumor samples
################################################################################
##################### COPYRIGHT ################################################
# New York Genome Center

# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright (2019) by the New York
# Genome Center. All rights are reserved. This software is supplied without
# any warranty or guaranteed support whatsoever. The New York Genome Center
# cannot be responsible for its use, misuse, or functionality.
# Version: 0.1
# Author: Kanika Arora (karora@nygenome.org) and Jennifer M Shelton
##################### /COPYRIGHT ###############################################
################################################################################
import sys
import os
import pandas as pd
##########################################################################
##############                  Custom functions              ############
##########################################################################


def prep_mantis(mantis, normal):
    '''
        Pre tables for mantis. Adjust threshold for
        projects without matched normals from
        0.4 to 0.62
        '''
    mantis_data = pd.read_csv(mantis, sep='\t')
    mantis_data = mantis_data[(mantis_data['Average Metric Value (Abbr)'] == 'Step-Wise Difference (DIF)')]
    mantis_table_data = mantis_data[['Status', 'Value', 'Threshold']].copy()
    if normal in ['NA12878']:
        mantis_table_data['Threshold'] = 0.62
        mantis_table_data.loc[(mantis_table_data.Value > mantis_table_data.Threshold), 'Status'] = 'Unstable'
        mantis_table_data.loc[(mantis_table_data.Value <= mantis_table_data.Threshold), 'Status'] = 'Stable'
    mantis_table_data.columns = ['MSI_Status', 'MSI_Score', 'Threshold']
    return mantis_table_data


def print_mantis(mantis_table_data, out):
    '''
        Write out Mantis CSV file
    '''
    mantis_table_data.to_csv(out, index=False, sep='\t')


def main():
    '''
        Update tables for mantis. Adjust threshold for
        projects without matched normals from
        0.4 to 0.62
    '''
    #  ==========================
    #  Input variables
    #  ==========================
    mantis = sys.argv[1]
    out = sys.argv[2]
    normal = sys.argv[3]
    assert os.path.isfile(mantis), 'Failed to find mantis file :' + mantis
    #  ==========================
    #  Run prep
    #  ==========================
    mantis_table_data = prep_mantis(mantis, normal)
    print_mantis(mantis_table_data, out)


##########################################################################
#####       Execute main unless script is simply imported     ############
#####                for individual functions                 ############
##########################################################################


if __name__ == '__main__':
    main()
