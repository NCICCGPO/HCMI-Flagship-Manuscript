#!/usr/bin/bash
set -e
cd pipelines/latent-tf-tumor-growth-distances
# Download data from GDC
bash scripts/gdc_download.sh
# Pre-process data
bash scripts/process.sh data/prep
# Calculate distances between sample pairs
python scripts/matrix_prep.py \
    --maindir data/prep \
    --outfile data/distance_metric/hcmi.counts.tsv
python scripts/matrix_intersect.py \
    -m data/distance_metric/hcmi.counts.tsv \
    -s src/PathwayCommons12.All.hgnc.sif.gz \
    --out data/distance_metric/features.txt
python scripts/matrix_expmax_normalize.py \
    data/distance_metric/hcmi.counts.tsv \
    --features data/distance_metric/features.txt \
    --out data/distance_metric/hcmi.normalized.tsv \
    --precision 5
python scripts/tf_net_train.py \
    data/distance_metric/hcmi.normalized.tsv src/PathwayCommons12.All.hgnc.sif.gz \
    --epochs 100 \
    -o data/distance_metric/model-tf-hcmi
python scripts/matrix_project.py \
    data/distance_metric/hcmi.normalized.tsv data/distance_metric/model-tf-hcmi \
    --out data/distance_metric/hcmi.latent.tsv
python scripts/matrix_distance.py \
    data/distance_metric/hcmi.latent.tsv \
    -m \
    -o data/distance_metric/hcmi.latent.dists.csv
python scripts/reformat_allvsall_lat.py \
    --dist data/distance_metric/hcmi.latent.dists.csv \
    --out data/distance_metric/main_results/latent_tf_dist_all_cohorts.tsv
cd ../..
