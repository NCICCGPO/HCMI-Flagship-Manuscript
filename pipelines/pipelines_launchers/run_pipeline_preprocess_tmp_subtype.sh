#!/usr/bin/bash
set -e
cancer=${1}
outdir=${2}
cd pipelines/classify-lab-models-and-tumors
# Download data from GDC
bash scripts/gdc_download.sh ${cancer}
# Pre-process data
bash scripts/process.sh ${cancer} ${outdir}
cd ../..