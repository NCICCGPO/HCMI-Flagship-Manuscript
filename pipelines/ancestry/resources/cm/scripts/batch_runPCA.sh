#!/bin/bash
script_root=`dirname $(readlink -f "$0")`
RSCRIPT=${script_root}/plotPCA.R

INPUTFILE=$1
ADMIXRESFILE=$2
OUTPUTPLOTFILE=$3

plink --bfile $INPUTFILE \
      --pca \
      --out $INPUTFILE

R --vanilla --args $INPUTFILE.eigenvec $ADMIXRESFILE $OUTPUTPLOTFILE <$RSCRIPT
