#!/usr/bin/env python
#	USAGE: python germ_filter_center_vcf.py VCF_FILE OUT_FILE
#   DESCRIPTION: filters out based on AF in Gnomad and removes svaba
################################################################################
##################### COPYRIGHT ################################################
# New York Genome Center

# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright (2023) by the New York
# Genome Center. All rights are reserved. This software is supplied without
# any warranty or guaranteed support whatsoever. The New York Genome Center
# cannot be responsible for its use, misuse, or functionality.
# Version: 0.2
# Author: Kanika Arora (karora@nygenome.org) and Jennifer M Shelton
##################### /COPYRIGHT ###############################################
################################################################################
import sys
import os
import logging as log
import pysam
from collections import Counter
##########################################################################
##############                  Custom functions              ############
##########################################################################

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


def read_vcf(vcf_file):
    '''
        Read in annotated VCF file.
    '''
    bcf_in = pysam.VariantFile(vcf_file)  # auto-detect input format
    return bcf_in

def get_csq_columns(bcf_in, source='CSQ'):
    '''
        get column names from the bar
        separated CSQ VEP annotation
        results. CSQ are Consequence
        annotations from Ensembl VEP.
    '''
    if source == 'CSQ':
        csq_columns = bcf_in.header.info[source].description.split()[-1].split('|') # grab the definitions
    else:
        csq_columns = bcf_in.header.info[source].description.split('Format:')[-1].split('|') # grab the definitions
        csq_columns = [col.strip() for col in csq_columns]
    return csq_columns


def get_csqs(record, csq_columns, source='CSQ'):
    '''
        Get new INFO field results.
    '''
    alt_count = len(record.alts)
    csq_dicts = {}
    spanning_deletion_offset = 0
    for i in range(alt_count):
        j = i - spanning_deletion_offset
        if record.alts[j] == '*':
            csq_dicts[i] = {csq_column : '' for csq_column in csq_columns}
            spanning_deletion_offset += 1
        else:
            try:
                csq_line = record.info[source][j]
            except UnicodeDecodeError: # for names with accents and other unexpected characters (rare)
                line = str(record)
                csq_line = line.split('\t')[7].split(source + '=')[1]
                csq_line = csq_line.split(';')[0]
                csq_line = csq_line.split(',')[j]
            csq_values = csq_line.split('|')
            csq_dict = dict(zip(csq_columns, csq_values))
            csq_dicts[i] = csq_dict
    return csq_dicts


def add_filter_header(bcf_out,
                      id,
                      description):
    '''
        Add new FILTER field
    '''
    bcf_out.header.filters.add(id=id,
                               number=None,
                               type=None,
                               description=description)
    return bcf_out

def get_af(af):
    '''
        return the greatest value recorded for AF.
    '''
    afs = [sing_af for sing_af in af.split('&') if not sing_af in ['', '.']]
    if len(afs) == 0:
        return ''
    try:
        return max([float(sing_af) for sing_af in afs])
    except ValueError:
        print([afs])
        print([af])
        sys.exit(1)

def test_af(record, alt_index, csq_columns, threshold=0.01):
    '''
        Test if AF > 0.01 in one germline database. PASS
        variants that don't have AF listed (for example
        in records from mouse databases)
        Pass variants that don't have sample columns (e.g. mouse variants from 00-All.normalized.vcf.gz)
    '''
    csq_dicts = get_csqs(record, csq_columns)
    GnomadGenomes_AF = get_af(csq_dicts[alt_index]['GnomadGenomes_AF_non_cancer'])
    GnomadExomes_AF = get_af(csq_dicts[alt_index]['GnomadExomes_non_cancer_AF'])
    AF_1000G = get_af(csq_dicts[alt_index]['AF_1000G'])
    afs = [db_af for db_af in [GnomadGenomes_AF, GnomadExomes_AF, AF_1000G] if not db_af == '']
    if any([float(db_af) > threshold for db_af in afs]):
        return True
    return False


def is_germline(record, alt, csq_columns, af=0.01):
    '''
        Check if matching variant is in the GRM VCF.
    '''
    for alt_index, alt in enumerate(record.alts):
        if test_af(record, alt_index, csq_columns, threshold=af):
            return True
    return False

def pull_svaba(record):
    ''''''
    # remove svaba
    if 'nygc_called_by' in record.info:
        nygc_called_by = record.info['nygc_called_by']
        if 'svaba' in nygc_called_by:
            # nygc record
            nygc_called_by = [c for c in nygc_called_by if not c == 'svaba']
            record.info['nygc_called_by'] = ','.join(nygc_called_by)
            nygc_num_callers = len(nygc_called_by)
            record.info['nygc_num_callers'] = nygc_num_callers
            # general record
            called_by = record.info['called_by']
            called_by = [c for c in called_by if not c == 'svaba']
            called_by = list((Counter(record.info['called_by']) - Counter(record.info['nygc_called_by'])).elements())
            record.info['called_by'] = called_by
            if nygc_num_callers < 2:
                return record, True
    return record, False

def pull_nygc_call(record):
    centers_called_by = [c for c in record.info['centers_called_by'] if not 'nygc' in c]
    record.info['centers_called_by'] = centers_called_by
    num_centers = len(centers_called_by)
    record.info['num_centers'] = num_centers
    if num_centers == 0:
        return record, True
    return record, False


def filter_vcf(bcf_in, bcf_out):
    '''
        Remove NYGC call if gnomad annotation would filter it out.
    '''
    germ_check = False
    if 'GnomadGenomes_AF_non_cancer' in bcf_in.header.info:
        germ_check = True
    svaba_check = False
    if 'nygc_called_by' in bcf_in.header.info:
        svaba_check = True
    csq_columns = get_csq_columns(bcf_in, source='CSQ')
    for record in bcf_in.fetch():
        remove_svaba = False
        if svaba_check:
            record, remove_svaba = pull_svaba(record)
        germ_filter = False
        skip_variant = False
        if germ_check:
            for alt in record.alts:
                if is_germline(record, alt, csq_columns):
                    germ_filter = True
        if germ_filter or remove_svaba:
            record, skip_variant = pull_nygc_call(record)
        if not skip_variant:
            # skip if only called by NYGC
            # and NYGC is removed
            exit_status = bcf_out.write(record)
            if exit_status != 0:
                print(exit_status)


def prep_header(bcf_in):
    if not 'num_centers' in bcf_in.header.info:
            bcf_in = add_info_header(bcf_out=bcf_in,
                              id='num_centers',
                              number='1',
                              type='Integer',
                              description='Number of centers that called the variant')
    if not 'centers_called_by' in bcf_in.header.info:
            bcf_in = add_info_header(bcf_out=bcf_in,
                              id='centers_called_by',
                              number='.',
                              type='String',
                              description='Name of the center(s) that the variant was called by')
    return bcf_in


def main():
    '''
        Change filter column from PASS to GRM or
        add GRM to filters.
    '''
    vcf_file = sys.argv[1]
    out_file = sys.argv[2]
    assert os.path.isfile(vcf_file), 'Failed to find somatic VCF call file :' + vcf_file
    bcf_in = read_vcf(vcf_file)
    bcf_in = prep_header(bcf_in)
    bcf_out = pysam.VariantFile(out_file, 'w', header=bcf_in.header)
    filter_vcf(bcf_in, bcf_out)


##########################################################################
#####       Execute main unless script is simply imported     ############
#####                for individual functions                 ############
##########################################################################
if __name__ == '__main__':
    main()
