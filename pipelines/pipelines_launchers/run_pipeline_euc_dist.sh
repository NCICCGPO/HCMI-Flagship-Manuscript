#!/usr/bin/bash
set -e
cancer=${1}
outdir=${2}
cd pipelines/euclidean-tumor-growth-distance 
# Download data from GDC
bash scripts/gdc_download.sh ${cancer}
# Pre-process data
bash scripts/process.sh ${cancer} ${outdir}
# Calculate distance between tumor and derived model
bash scripts/euc_distance.sh ${cancer} ${outdir}
cd ../..
