#!/bin/bash

# script to call admixture program
# takes as input:
# file path/root for PLINK-formatted genotypes
#    assumes this includes .bed,.bim,.fam
#    assumes in a folder that one is allowed to write to
# FOR NOW:
# assumes genotypes were merged with 1kg reference samples
# assumes running to estimate proportion ancestry in 5 superpopulations only

while [[ $# -gt 1 ]]; do
    key=${1}
    case ${key} in
    -script_root)
        script_root=${2};;
    -plink)
        plink_input=${2};;
    -super_pops)
        super_pops=${2};;
    *)
        echo "Unknown option: ${key}."
        echo "Usage: -plink <plink> -super_pops <super population>"
        exit 65
    ;;
    esac
    shift
    shift
done

RSCRIPT_INPUT=${script_root}/../resources/cm/scripts/formatADMIXTUREinputs.subpops.R
RSCRIPT_OUTPUT=${script_root}/../resources/cm/scripts/formatADMIXTUREoutputs.R
kgpops=${script_root}/../resources/cm/1000Genomes/1kgPops.txt
super_pops=${super_pops:-superpops}

INPUT=`readlink -f ${plink_input}`

R --vanilla --args $INPUT $kgpops $super_pops <${RSCRIPT_INPUT}

cd $(dirname $INPUT)

numpops=`cut -d " " -f 1 ${INPUT}.fam | sort | uniq | grep -v 0 | wc -l`
echo admixture --supervised $INPUT.bed $numpops

admixture --supervised $INPUT.bed $numpops

paste $INPUT.pop \
      $INPUT.fam \
      $INPUT.$numpops.Q \
      > $INPUT.admixtureRes.subpops

# now, read in file to R and remove all non-study sample``
# also name the columns appropriately
R --vanilla --args $INPUT.admixtureRes.subpops <${RSCRIPT_OUTPUT}
