#!/usr/bin/env python2.7

# New York Genome Center
# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright (2014) by the New York
# Genome Center. All rights are reserved. This software is supplied without
# any warranty or guaranteed support whatsoever. The New York Genome Center
# cannot be responsible for its use, misuse, or functionality.
# Version: 0.1
# Author: Kanika Arora (karora@nygenome.org)
# Adapted from verion 0.4 of calculate_merge_project_genotypes.py 
# Author: Ewa A Grabowska (egrabowska@nygenome.org)

import sys
import os
import optparse
import math
import imp
import glob
from collections import defaultdict
import json
from tempfile import NamedTemporaryFile
from shutil import move

desc = """Program to generate genotype matrix for a sample based on pileup file"""
parser = optparse.OptionParser(version='%prog version 0.5 16/January/2019', description=desc)
parser.add_option('-P', '--pileup_file', help='Input pileup file [required]', action='store')
parser.add_option('-G', '--genomebuild', help='Reference genome build', action="store",default='b38',choices=['b37','b38'])
parser.add_option('-T', '--type', help='Choose whether the sample data is whole genome (wgs) or exome sequencing (exome)', default="wgs",choices=['wgs','exome','impact','rna'])
parser.add_option('-C', '--min_cov', help='MIN COVERAGE TO CALL GENOTYPE', default=10, type='int', action='store')
parser.add_option('-O', '--outfile', help='TXT OUTPUT FILE [required]', type='string', action='store')
parser.add_option('-Q', '--min_mapping_quality', help='MIN MAPPING QUALITY', default=10, type='int', action='store')
parser.add_option('-B', '--min_base_quality', help='MIN BASE QUALITY', default=20, type='int', action='store')
parser.add_option('-R', '--repository', help='Directory with required modules scripts', default=os.path.dirname(os.path.abspath(sys.argv[0])))

(opts, args) = parser.parse_args()
pathname = os.path.dirname(sys.argv[0])
markers_json_file=os.path.join(pathname,"../markers/path_to_reference_and_markers_for_genotyping.json")
json_data=open(markers_json_file).read()
data = json.loads(json_data)

if opts.repository:
    DIR=opts.repository

sys.path.append(DIR)

ContaminationMarker = imp.load_source('/ContaminationMarker', DIR + '/ContaminationMarker.py')

if not opts.pileup_file or not opts.outfile:
    parser.print_help()
    sys.exit(1)

if opts.genomebuild:
    BUILD = opts.genomebuild

if opts.type:
    TYPE = opts.type

MARKER_FILE = data[BUILD][TYPE]["markerFilePrefix"]+".markers.txt"
REFERENCE = data[BUILD]["reference"]
REMOVE_CHR_PREFIX = data[BUILD]["remove_chr_prefix"]

if not os.path.exists(MARKER_FILE):
    print('ERROR: Marker file {0} cannot be found.'.format(MARKER_FILE))
    sys.exit(2)
    
if not os.path.exists(opts.pileup_file):
    print('ERROR: Input pileup file {0} cannot be found.'.format(opts.pileup_file))
    sys.exit(2)
    
outfile = open(opts.outfile, 'w')
    
Markers = ContaminationMarker.get_markers(MARKER_FILE)
COVERAGE_THRESHOLD = opts.min_cov
MMQ = opts.min_mapping_quality
MBQ = opts.min_base_quality
##AA_BB_only = opts.normal_homozygous_markers_only
        
genotype_likelihoods = ContaminationMarker.genotype_likelihoods_for_markers(Markers, opts.pileup_file, min_map_quality=MMQ, min_base_quality=MBQ)
D = genotype_likelihoods


G = ['0/0', '0/1', '1/1']

for marker in Markers:
    line_to_output = "\t".join([Markers[marker].chrom, Markers[marker].pos, Markers[marker].id, Markers[marker].ref, Markers[marker].alt])
    GL = D[marker]
    if GL is None or GL['coverage'] < COVERAGE_THRESHOLD:
        line_to_output += "\tNA"
	#continue
    else:
        genotype = G[GL['likelihoods'].index(max(GL['likelihoods']))]
        line_to_output += "\t" + genotype
    outfile.write(line_to_output + "\n")
outfile.close()

#  ===============================
#  Update chromosome prefix if
#  requested
#  ===============================
if REMOVE_CHR_PREFIX == "yes":
    print("Removing 'chr' prefix...")

    with NamedTemporaryFile(delete=False) as tmp_source:
        with open(opts.outfile) as source_file:
            for line in source_file:
                if line.startswith("chr"):
                    tmp_source.write(line[3:])
                else:
                    tmp_source.write(line)

    move(tmp_source.name, source_file.name)

#  ===============================
#  Set output permissions
#  ===============================
try:
    os.chmod(opts.outfile, 0660)
except OSError:
    pass # only the owner can change the file permissions
