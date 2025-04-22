#!/usr/bin/bash
set -e
cancer=${1}
outdir=${2}
datatype=${3}
cd pipelines/classify-lab-models-and-tumors
# Run Gene Expression Classifier
if [[ ${datatype} == 'GEXP' ]]; then
    # where specify cancer, tumor-file, model-file, transformed-dir
    bash scripts/run_classify_GEXP.sh \
        ${cancer} \
        ${outdir}/${cancer}_GEXP/${cancer}_GEXP_prep2_Tumor.tsv \
        ${outdir}/${cancer}_GEXP/${cancer}_GEXP_prep2_Model.tsv \
        data/classifier_gexp/ml_ready_qrank
elif [[ ${datatype} == 'METH' ]]; then
    # where specify cancer, tumor-file, model-file, transformed-dir
    bash scripts/run_classify_METHYL.sh \
        ${cancer} \
        data/classifier_methyl/processed/HCMI_TMP_subtype_prediction_feature_matrix_${cancer}.tsv
fi
cd ../..
