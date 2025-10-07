#!/bin/bash
set -euo pipefail
SRCDIR=$(realpath $(dirname $0)) 



## Dirs
PROJECT_DIR=/path/to/directory/above/output/dir                     ## User choice 
JABBA_DIR=$PROJECT_DIR/path/to/jabba/output                         ## User will need to obtain
COVERAGE_DIR=$PROJECT_DIR/path/to/fragcounter/output                ## User will need to obtain

AA_DIR=$PROJECT_DIR/amplicon_architect
WORKING_DIR=$AA_DIR/tumor-model-comparisons
OUT_DIR=$WORKING_DIR/feature-similarity
FEAT_DIR=$OUT_DIR/features-to-graph
FIG_DIR=$WORKING_DIR/fig
LOG_DIR=$WORKING_DIR/logs


## Metadata
METADATA_DIR=$SRCDIR/../metadata                                  
TUMOR_MODEL_PAIRS=$METADATA_DIR/tumor_model_pairs_all.txt           
TNFILE_CSV=$METADATA_DIR/older_pairings/official_pairs.csv          
BAM_MAP=$METADATA_DIR/nygc_bam_map.csv                              
PURITY_PLOIDY=$METADATA_DIR/2.15.24_consensus_pp.txt

## Software
AA_REPO=/path/to/AmpliconClassifier                                 ## User needs to fill this out 


## Parameters
REFERENCE_BUILD=GRCh38
SIMILARITY_MINIMUM=0.75



## Make output dirs
mkdir -p \
    $OUT_DIR \
    $FEAT_DIR \
    $FIG_DIR \
    $LOG_DIR \
    $AA_DIR/bed/tmp \
    $AA_DIR/input-tables \
    $AA_DIR/logs \
    $AA_DIR/joint-calling \
    $AA_DIR/reconstructions




############################### 
## Post-process JaBbA output ##
###############################

## Export from VCF
while read participantId tumor_aliquot_barcode normal_aliquot_barcode discovery pairId tumorBam normalBam group_normal_aliquot_barcode group_tumor_aliquot_barcode linkage_file_type; do
    
    tn=$tumor_aliquot_barcode--$normal_aliquot_barcode

    in_file=$JABBA_DIR/$tn/jabba.simple.cnv.vcf
    out_file=$AA_DIR/bed/tmp/$tn.jabba.cnv.bed

    if [[ ! -f $in_file ]]; then 
        >&2 echo "Missing jabba output for $tn"
        continue 
    fi

    cmd="Rscript $SRCDIR/cnv-vcf-to-bed.r -i $in_file -o $out_file"
    jobname=jabba_bed.$tn

    sbatch \
        -o $LOG_DIR/$jobname.o.%A \
        -c 1 \
        --mem=4G \
        --job-name=$jobname \
        --wrap="$cmd"

done < <(tail -n +2 $TNFILE_CSV | tr ',' '\t')


## Trim BED to how AA expects 
while read participantId tumor_aliquot_barcode normal_aliquot_barcode discovery pairId tumorBam normalBam group_normal_aliquot_barcode group_tumor_aliquot_barcode linkage_file_type; do
    
    tn=$tumor_aliquot_barcode--$normal_aliquot_barcode

    in_file=$AA_DIR/bed/tmp/$tn.jabba.cnv.bed
    out_file=$AA_DIR/bed/$tn.jabba.bed

    if [[ ! -f $in_file ]]; then 
        >&2 echo "Missing jabba output for $tn"
        continue 
    fi

    tail -n +2 $in_file | cut -f 1-4 | xargs -I {} echo "chr{}"> $out_file

done < <(tail -n +2 $TNFILE_CSV | tr ',' '\t')




####################################################
## Run Amplicon Architect in "joint calling" mode ##
####################################################

## Build input tables
Rscript $SRCDIR/init-joint-calling-input-tables.r \
    --tn_file $TNFILE_CSV \
    --cnv_dir $AA_DIR/bed \
    --bam_map $BAM_MAP \
    --out_dir $AA_DIR/input-tables

CMD_FILE=$AA_DIR/task_array_cmds.sh


## Build task array 
while read participantId; do 

    in_file=$AA_DIR/input-tables/$participantId.joint_aa_input.txt
    out_dir=$AA_DIR/joint-calling/$participantId

    if [[ ! -f $in_file ]]; then 
        >&2 echo "Missing input file for $participantId"
        continue 
    fi 

    mkdir -p $out_dir
    rm -r $out_dir/* || true 

    cmd="/nfs/sw/amplicon/amplicon-1.15.2/GroupedAnalysisAmpSuite.py \
            -i $in_file \
            -o $out_dir \
            -t 4 \
            --ref $REFERENCE_BUILD"

    echo $cmd

done < <(tail -n +2 $TNFILE_CSV | cut -f1 -d, | sort -u) > $CMD_FILE


## Kick off task array 
cmd="bash $SRCDIR/run-job-array-task.sh $CMD_FILE"
jobname=hcmi_joint_aa_pipeline
njobs=$(cat $CMD_FILE | wc -l)

sbatch \
    -o $AA_DIR/logs/$jobname.o.%A.%a \
    --array [1-$njobs]%100 \
    -c 4 \
    --mem=20G \
    --job-name=$jobname \
    --wrap="$cmd"



## Overall AA summary
Rscript $SRCDIR/init-summarize-aa-results.r \
    --tn_file $TNFILE_CSV \
    --tm_file $TUMOR_MODEL_PAIRS \
    --aa_dir $AA_DIR/joint-calling \
    --out_file $FIG_DIR/aa-amplicon-size-summary.pdf \
    --out_file_txt $WORKING_DIR/aa-amplicon-summary-merged.txt \
    --out_file_txt_unfiltered $WORKING_DIR/aa-amplicon-summary-merged.unfiltered.txt




########################
## Feature comparison ##
########################

## For each tumor model pair
while read model tumor; do 

    ## Pull normal from tumor normal pair list
    normal=$(awk -F, -v t=$tumor '($2 == t) {print $3}' $TNFILE_CSV)
    patient=$(awk -F, -v t=$tumor '($2 == t) {print $1}' $TNFILE_CSV)

    if [[ -z $normal ]]; then
        continue
    fi

    tn=$tumor--$normal
    mn=$model--$normal

    ## Find feature files if they exist and concatenate
    tn_in_file=$AA_DIR/joint-calling/$patient/$tn/${tn}_classification/${tn}_features_to_graph.txt
    mn_in_file=$AA_DIR/joint-calling/$patient/$mn/${mn}_classification/${mn}_features_to_graph.txt

    tm_feature_file=$FEAT_DIR/${tumor}--${model}_features_to_graph.txt

    rm $tm_feature_file || true
    [[ -f $tn_in_file ]] && cat $tn_in_file >> $tm_feature_file
    [[ -f $mn_in_file ]] && cat $mn_in_file >> $tm_feature_file


    ## If we have a feature file, run comparison
    if [[ -f $tm_feature_file ]]; then

        cmd="$AA_REPO/feature_similarity.py \
            -f $tm_feature_file \
            -o $OUT_DIR/$tumor--$model \
            --ref GRCh38"

        jobname=featsim.$tumor--$model

        sbatch \
            -o $LOG_DIR/$jobname.o.%A \
            --mem=4G \
            -c 1 \
            --job-name=$jobname \
            --wrap="$cmd"

    fi

done < <(tail -n +2 $TUMOR_MODEL_PAIRS)


## Merge and annotate amplicon comparisons
Rscript $SRCDIR/init-merge-and-annotate-feature-comparisons.r \
    --tm_file $TUMOR_MODEL_PAIRS \
    --aa_dir $AA_DIR/joint-calling \
    --tn_file $TNFILE_CSV \
    --fc_dir $OUT_DIR \
    --min_jaccard_similarity $SIMILARITY_MINIMUM \
    --amplicons $WORKING_DIR/aa-amplicon-summary-merged.unfiltered.txt \
    --out_file $WORKING_DIR/hcmi-merged-feature-comparisons.txt


## Rescue filtered ecDNAs if they were found in paired tumor/normal 
Rscript $SRCDIR/init-rescue-ecdna-calls.r \
    --in_file_amplicons $WORKING_DIR/aa-amplicon-summary-merged.unfiltered.txt \
    --in_file_concordant_ecdna $WORKING_DIR/hcmi-merged-feature-comparisons.concordant_ecdna.txt \
    --out_file $WORKING_DIR/aa-amplicon-summary-merged.rescued_ecdna.txt




###############################################
## Try to generate candidate reconstructions ##
###############################################

while read pair_id amp_id; do 

    sample_id=$(echo $pair_id | sed -e 's/--.*//g')
    patient=$(awk -F, -v t=$sample_id '($2 == t) {print $1}' $TNFILE_CSV)

    feat2graph=$AA_DIR/joint-calling/$patient/$pair_id/${pair_id}_classification/${pair_id}_features_to_graph.txt
    class_bed=$AA_DIR/joint-calling/$patient/$pair_id/${pair_id}_classification/${pair_id}_classification_bed_files/${amp_id}_intervals.bed

    amp_graph=$(awk -v b=$class_bed '$1 == b {print $2}' $feat2graph)

    cmd="/nfs/sw/amplicon/amplicon-1.15.2/scripts/CAMPER.py \
        -g $amp_graph \
        --remove_short_jumps \
        --runmode bulk \
        --outname $AA_DIR/reconstructions/$amp_id.reconstruction"

    jobname=$amp_id.ecdna_reconstruction

    sbatch \
        -o $LOG_DIR/$jobname.o.%A \
        --mem=8G \
        -c 1 \
        --job-name=$jobname \
        --wrap="$cmd"

done < <(tail -n +2 $WORKING_DIR/aa-amplicon-summary-merged.rescued_ecdna.txt | awk -F "\t" '$4 == "ecDNA" {print $1, $3}')




#########################
## Extended Data plots ##
#########################

## Extended Data 3d: concordance barplot split by binned purity
Rscript $SRCDIR/plotting/plt-extended-data-3c-concordance-by-purity.r \
    --tm_file $TUMOR_MODEL_PAIRS \
    --tn_file $TNFILE_CSV \
    --amplicons $WORKING_DIR/aa-amplicon-summary-merged.rescued_ecdna.txt \
    --purity_ploidy $PURITY_PLOIDY \
    --out_file $FIG_DIR/extended-data-3c-concordance-by-purity.svg


## Extended Data 9e: MYCN ecDNA coverage gTrack for HCM-BROD-0613-C71
## NOTE: requires fragcounter coverage 
Rscript $SRCDIR/plotting/plt-extended-data-9e-paired-coverage-tracks.r \
            --tumor HCM-BROD-0613-C71-01A \
            --model HCM-BROD-0613-C71-85A \
            --tumor_coverage $COVERAGE_DIR/$tumor/cov.rds \
            --model_coverage $COVERAGE_DIR/$model/cov.rds \
            --amp1 HCM-BROD-0613-C71-85A-01D-A80U-36--HCM-BROD-0613-C71-10A-01D-A80U-36_amplicon1_ecDNA_2 \
            --amp2 HCM-BROD-0613-C71-01A-11D-A80U-36--HCM-BROD-0613-C71-10A-01D-A80U-36_amplicon1_ecDNA_1 \
            --field reads.corrected \
            --feat2graph=$tm_feature_file \
            --genes MYCN \
            --out_file=$FIG_DIR/extended-data-3c-mycn-amp.svg
