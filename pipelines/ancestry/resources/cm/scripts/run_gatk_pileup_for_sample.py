#!/usr/bin/env python2.7

# New York Genome Center
# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright (2014) by the New York
# Genome Center. All rights are reserved. This software is supplied without
# any warranty or guaranteed support whatsoever. The New York Genome Center
# cannot be responsible for its use, misuse, or functionality.
# Version: 0.3
# Author: Ewa A Grabowska (egrabowska@nygenome.org)
# Minor changes made by
# Kanika Arora (karora@nygenome.org)
# 2018-03-14
import os
import sys
import optparse
from shutil import move
from tempfile import NamedTemporaryFile
import json

desc = """Program to run GATK Pileup on a single sample"""
parser = optparse.OptionParser(version='%prog version 0.3 14/June/2016', description=desc)
parser.add_option('-B', '--bam', help='BAMFILE [mandatory field]', action='store')
parser.add_option('-O', '--outfile', help='OUTPUT FILE (PILEUP) [mandatory field]', type='string', action='store')
parser.add_option('-G', '--genomebuild', help='Reference genome build', action="store",default='b38',choices=['b37','b38'])
parser.add_option('-T', '--type', help='Choose whether the sample data is whole genome (wgs) or exome sequencing (exome)', default="wgs",choices=['wgs','exome',"rna","impact"])
parser.add_option('-J', '--gatk', help='GATK JAR', action='store')
#parser.add_option('--remove_chr_prefix', help='REMOVE CHR PREFIX FROM THE CHROMOSOME COLUMN IN THE OUTPUT FILE [false by default]', default=False, action='store_true')

pathname = os.path.dirname(sys.argv[0])
markers_json_file=os.path.join(pathname,"../markers/path_to_reference_and_markers_for_genotyping.json")
json_data=open(markers_json_file).read()
data = json.loads(json_data)


(opts, args) = parser.parse_args()

if not opts.bam or not opts.outfile:
    parser.print_help()
    sys.exit(1)
    
if not os.path.exists(opts.bam):
    print('ERROR: Specified bamfile {0} cannot be found.'.format(opts.bam))
    sys.exit(1)
    
GATK = opts.gatk

if not os.path.exists(GATK):
    print('ERROR: GATK jar {0} cannot be find.'.format(GATK))
    sys.exit(2)   

if opts.genomebuild:
    BUILD = opts.genomebuild

if opts.type:
    TYPE = opts.type

MARKER_FILE = data[BUILD][TYPE]["markerFilePrefix"]+".markers.bed"
REFERENCE = data[BUILD]["reference"]

if not os.path.exists(MARKER_FILE):
    print('ERROR: Marker file {0} cannot be found.'.format(MARKER_FILE))
    sys.exit(2)
    
if not os.path.exists(REFERENCE):
    print('ERROR: Reference genome {0} cannot be found.'.format(REFERENCE))
    sys.exit(3)
#  ===============================
#  Call GATK
#  ===============================
command_line = "java -Xmx12g -jar {0} -T Pileup -R {1} -I {2} -L {3} -o {4} -verbose -rf DuplicateRead --filter_reads_with_N_cigar --filter_mismatching_base_and_quals".format(GATK, REFERENCE, opts.bam, MARKER_FILE, opts.outfile)
print command_line
os.system(command_line)

#  ===============================
#  Set output permissions
#  ===============================
try:
    os.chmod(opts.outfile, 0660)
except OSError:
    pass # only the owner can change the file permissions
