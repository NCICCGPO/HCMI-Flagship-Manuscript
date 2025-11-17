#!/usr/bin/env python

"""
Script to generate *.toGenos file from GVCF.
"""

import os
import argparse
import logging

from cyvcf2 import VCF

# Setup Argparse
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('-i', dest='gvcf', required=True)
parser.add_argument('-m', dest='markers', required=True)
parser.add_argument('-d', dest='dbsnp', required=True)
parser.add_argument('-r', dest='ref_ver', required=True, choices=['b37', 'b38'])
parser.add_argument('-o', dest='output', required=True)
parser.add_argument('-v', dest='verbose', action='store_true')
args = parser.parse_args()

# Setup logger
log = logging.getLogger(os.path.basename(__file__))
if args.verbose:
    log.setLevel(logging.DEBUG)
else:
    log.setLevel(logging.INFO)
ch = logging.StreamHandler()
formatter = logging.Formatter(
    '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )
ch.setFormatter(formatter)
log.addHandler(ch)


class VcfRecord():
    """
    class to represent a record in the vcf file.
    """

    def __init__(self, cyvcf2_var_rec):
        self.parse_cyvcf2_record(cyvcf2_var_rec)

    def parse_cyvcf2_record(self, rec):
        self.chrom = rec.CHROM
        self.pos = rec.end
        self.id = rec.ID
        self.ref = rec.REF
        self.alt = rec.ALT
        self._rec = rec

    def __repr__(self):
        return self._rec

    def isIndel(self):
        # Deletion
        if len(self.ref) > 1:
            return True
        for a in self.alt:
            # Ignore reference symbolic allele
            if a == '<NON_REF>' or a == []:
                continue
            # Insertion
            elif len(a) > 1:
                return True
        return False

    def isRef(self):
        if len(self.alt) == 1 and '<NON_REF>' in self.alt:
            return True
        elif self.alt == []:
            return True
        return False

    def isPosMatch(self, pos):
        if self.pos == pos:
            return True
        return False

    def isConcordant(self, vcf_record):
        if self.ref == vcf_record.ref and self.alt[0] in vcf_record.alt:
            return True
        return False

    def getGT(self):
        # HET
        if self._rec.gt_types[0] == 1:
            return '1'
        # HOM_REF
        elif self._rec.gt_types[0] == 0:
            return '2'
        # HOM_ALT
        elif self._rec.gt_types[0] == 3:
            return '0'
        else:
            # print self._rec
            raise RuntimeError('Not expecting GT: %s' % (self._rec.gt_types))

def query_gvcf(vcf, chrom, pos, log):
    """
    Query VCF and return record if coordinate position
    matches and whether it is ref block, else return None.
    return: is_gvcf_block, VcfRecord
    """
    for v in vcf('%s:%s-%s' % (chrom, pos, pos)):
        log.debug(v)
        vr = VcfRecord(v)
        try:
            vr.getGT()
        except:
            continue
        if vr.isIndel():
            continue
        elif vr.isRef():
            return (True, vr)
        elif vr.isPosMatch(pos):
            return (False, vr)
    return (None, None)

def query_vcf(vcf, chrom, pos, log):
    """
    Query VCF and return record if coordinate position
    matches and whether it is ref block, else return None.
    return: is_gvcf_block, VcfRecord
    """
    for v in vcf('%s:%s-%s' % (chrom, pos, pos)):
        log.debug(v)
        vr = VcfRecord(v)
    #    try:
    #        vcf_rec.getGT()
    #    except:
    #        continue
        if vr.isIndel():
            continue
        elif vr.isRef():
            return (True, vr)
        elif vr.isPosMatch(pos):
            return (False, vr)
    return (None, None)


def main(args):
    log.info("Reference version: %s" % args.ref_ver)

    gvcf = VCF(args.gvcf)
    dbsnp = VCF(args.dbsnp)

    sample_name = gvcf.samples[0]
    sample_name=sample_name.replace("_","-")
    HEADER = "CHR\tSNP\t(C)M\tPOS\tCOUNTED\tALT\t0_%s" % sample_name
    outfh = open(args.output, 'w')
    outfh.write(HEADER + '\n')

    # Total markers
    t_cnt = 0
    # Markers skipped
    m_cnt = 0

    with open(args.markers, 'r') as fh:
        for line in fh.readlines():
            log.debug('Marker: %s' % line)
            t_cnt += 1
            row = line.strip().split()
            chrom = row[0]
            pos = int(row[1])
            rsid = row[2]
            output = [
                chrom,
                rsid,
                '0',
                str(pos)
            ]

            vcf_ref_block, vcf_rec = query_gvcf(gvcf, chrom, pos, log)
            if args.ref_ver == 'b38':
                dbsnp_ref_block, dbsnp_rec = query_vcf(dbsnp, chrom.lstrip('chr'), pos, log)
            else:
                dbsnp_ref_block, dbsnp_rec = query_vcf(dbsnp, chrom, pos, log)

            # Unable to find the position in dbSNP vcf
            if dbsnp_rec is None:
                log.info('dbSNP record not found, skipping: %s' % rsid)
                m_cnt += 1
                continue
            else:

                a1 = dbsnp_rec.ref
                a2 = dbsnp_rec.alt[0]

                if a1 == a2:
                    msg = "dbSNP file %s a1 == a2 %s" % args.dbsnp
                    log.error(msg)
                    raise RuntimeError(msg)

            # No entry in VCF
            if vcf_rec is None:
                log.debug('VCF record not found: %s' % rsid)
                output.extend([a1, a2, 'NA'])

            # Record is reference block
            elif vcf_ref_block:
                output.extend([a1, a2, '2'])

            # check if ref and alt bases match
            elif vcf_rec.isConcordant(dbsnp_rec):
                output.extend([a1, a2, vcf_rec.getGT()])

            # Coordinate position matches but alt allele does not match
            else:
                output.extend([a1, a2, 'NA'])

            outfh.write('\t'.join(output) + '\n')

        log.info("Total markers: %d" % t_cnt)
        log.info("Markers skipped: %d" % m_cnt)
        outfh.close()

main(args)
