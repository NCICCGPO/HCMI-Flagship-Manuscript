#!/usr/bin/env python


import argparse
import pysam


def modify_header(bcf_in):
    '''
        Add new FORMAT fields
    '''
    bcf_in.header.formats.add(id='GT',
                           number='1',
                           type='String',
                           description='Genotype.')
    return bcf_in


def read_vcf(vcf_file):
    '''
        Read in annotated VCF file.
    '''
    bcf_in = pysam.VariantFile(vcf_file)  # auto-detect input format
    return bcf_in


def get_args():
    parser = argparse.ArgumentParser(
                                     description='DESCRIPTION: Add GT if missing.')
                                     #    Documentation parameters
    parser.add_argument('-v', '--vcf',
                        dest='vcf',
                        help='Input VCF file',
                        required=True)
    parser.add_argument('-o', '--out',
                        dest='out',
                        help='Output VCF file',
                        required=True)
    parser.add_argument('-n', '--normals',
                        dest='normals',
                        nargs='*',
                        help='Normal sample ids',
                        required=True)
    args = parser.parse_args()
    return args.__dict__

def add_record_gt(record, samples, normals):
    formats = record.format.keys()
    if not 'GT' in formats:
        for sample in samples:
            if sample in normals:
                record.samples[sample]['GT'] = (0,0)
            else:
                record.samples[sample]['GT'] = (0,1)
    return record

def add_gt():
    '''
    Add GT if missing
    '''
    args = get_args()
    normals = args['normals']
    bcf_in = read_vcf(args['vcf'])
    samples = list(bcf_in.header.samples)
    formats = bcf_in.header.formats.keys()
    if 'GT' not in formats:
        bcf_in = modify_header(bcf_in)
    bcf_out = pysam.VariantFile(args['out'], 'w', header=bcf_in.header)
    for record in bcf_in.fetch():
        if 'GT' not in formats:
            record = add_record_gt(record, samples, normals)
        exit_status = bcf_out.write(record)
        if exit_status != 0:
            print(exit_status)



def main():
    add_gt()

##########################################################################
#####       Execute main unless script is simply imported     ############
#####                for individual functions                 ############
##########################################################################


if __name__ == '__main__':
    main()