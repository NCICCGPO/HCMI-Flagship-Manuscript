#!/bin/bash
set -euxo pipefail
source /etc/profile.d/modules.sh
SRCDIR=$(realpath $(dirname $0))
module purge 
module load amplicon/1.15.2 \
            bedtools/2.31.1

AA_THREADS=4
AA_DATA_REPO=/gpfs/commons/projects/TCGA/gdc-awg/HCMI/WGS/compbio/amplicon_architect/data_repo
AC_SRC=/nfs/sw/amplicon/amplicon-1.15.2/AmpliconClassifier
export AC_SRC AA_DATA_REPO

## Use personal mosek license
MOSEKLM_LICENSE_FILE=/gpfs/commons/home/whooper/mosek
export MOSEKLM_LICENSE_FILE

## Run command 
CMD_LIST=$1


## Run for this sample
while read cmd; do

    eval "$cmd"

done < <(head -n $SLURM_ARRAY_TASK_ID $CMD_LIST | tail -n 1)