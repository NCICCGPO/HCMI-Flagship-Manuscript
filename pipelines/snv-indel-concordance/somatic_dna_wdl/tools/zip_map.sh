#!/bin/bash

# USAGE: bash zip_map.sh BamMapLike

# DESCRIPTION: zip two named arrays and output the list items as key
# value pairs

BamMapLike=$1

set -e -o pipefail

cat ${BamMapLike} | jq '[.sampleId, .bam] | transpose | map( {(.[0]) :  .[1]} ) | add '