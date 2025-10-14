#!/usr/bin/env python
#	USAGE: python filter_baf.py VCF OUT
#   DESCRIPTION: Filter VCF for lines with a rsID in the Existing variation annotation from VEP
################################################################################
##################### COPYRIGHT ################################################
# New York Genome Center

# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright (2019) by the New York
# Genome Center. All rights are reserved. This software is supplied without
# any warranty or guaranteed support whatsoever. The New York Genome Center
# cannot be responsible for its use, misuse, or functionality.
# Author: Kanika Arora (karora@nygenome.org) and Jennifer M Shelton
##################### /COPYRIGHT ###############################################
################################################################################
import sys
import os
import pysam
##########################################################################
##############                  Custom functions              ############
##########################################################################


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


def get_dbsnp_rs(Existing_variation):
    '''
        Remove Cosmic IDs and split by comma.
    '''
    ids = Existing_variation.split('&')
    good_ids = [id for id in ids if id.startswith('rs')]
    return good_ids


def get_rsID(csq_dicts):
    '''
        Get dbSNP annotation.
    '''
    dbSNP_RSs = get_dbsnp_rs(csq_dicts['Existing_variation'])
    return dbSNP_RSs


def test_af(af, threshold=0.65):
    '''
    test that the AF is greater than X%
    '''
    if af == '':
        return False
    elif float(af) > threshold:
        return True
    return False


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


def filter_novel(record, csq_columns, cgc_columns, threshold=0.65):
    '''
        Filter rows that do not include
        variants found in Gnomad with a AF > 65%.
        Lines are also kept if they have:
            'HIGH' or 'MODERATE' impact and
            in CGC (either tier 1 or 2 is fine)
    '''
    csq_dicts = get_csqs(record, csq_columns)
    if cgc_columns:
        cgc_dicts = get_csqs(record, cgc_columns, source='CancerGeneCensus')
    else:
        cgc_dicts = False
    for i, alt in enumerate(record.alts):
#        dbSNP_RSs = get_rsID(csq_dicts[i])
        AF_1000G = get_af(csq_dicts[i]['AF_1000G'])
        GnomadGenomes_AF = get_af(csq_dicts[i]['GnomadGenomes_AF'])
        GnomadExomes_AF = get_af(csq_dicts[i]['GnomadExomes_AF'])
        if test_af(GnomadGenomes_AF, threshold=threshold) \
            or test_af(GnomadExomes_AF, threshold=threshold):
            return True
        impact = csq_dicts[i]['IMPACT']
        if cgc_dicts:
            tier = cgc_dicts[i]['Tier']
            if impact in ['HIGH', 'MODERATE'] and \
                    tier in ['1', '2']:
                return True
    return False


def write_file(bcf_in, vcf_out_file, csq_columns, cgc_columns,
               threshold=0.65):
    '''
        Write out the header
    '''
    #  =====================
    #  Write out split VCF
    #  =====================
    bcf_out = pysam.VariantFile(vcf_out_file, 'w', header=bcf_in.header)
    for record in bcf_in.fetch():
        if filter_novel(record, csq_columns,
                        cgc_columns, threshold=threshold):
            exit_status = bcf_out.write(record)
            if exit_status != 0:
                print(exit_status)
    return True


def main():
    '''
        Filter VCF for lines with a dbSNP rsID
    '''
    #  ==========================
    #  Input variables
    #  ==========================
    vcf_file = sys.argv[1]
    vcf_out_file = sys.argv[2]
    if len(sys.argv) > 3:
        threshold = float(sys.argv[3])
    else:
        threshold = 0.65
    assert os.path.isfile(vcf_file), 'Failed to find VCF file :' + vcf_file
    #  ==========================
    #  Run filter
    #  ==========================
    bcf_in = read_vcf(vcf_file)
    csq_columns = get_csq_columns(bcf_in)
    if 'CancerGeneCensus' in bcf_in.header.info:
        cgc_columns = get_csq_columns(bcf_in, source='CancerGeneCensus')
    else:
        cgc_columns = False
    write_file(bcf_in, vcf_out_file, csq_columns,
               cgc_columns, threshold=threshold)



##########################################################################
#####       Execute main unless script is simply imported     ############
#####                for individual functions                 ############
##########################################################################


if __name__ == '__main__':
    main()
