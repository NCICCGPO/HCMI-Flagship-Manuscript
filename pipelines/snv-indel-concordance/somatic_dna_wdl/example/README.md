# Tutorial



## WGS Human Somatic Pipeline v7
- [General notes](#notes)
- [Run full pipeline from BAM](#from_bam)
- [Run SNV, INDEL, SV and CNV calling from BAM](#snv_indel_from_bam)
- [Run full pipeline from FASTQ](#from_fastq)


### General notes
<a name="notes"></a>

BAI files should end in the extension `.bai`

Bams can either be passed in the input JSON as fully formed objects (see below):

```
"pairInfos": [
        {
            "tumor": "TUMOR",
            "normal": "NORMAL",
            "normalFinalBam": {
                "bam": " gs://YOUR_BUCKET/NORMAL.bam",
                "bamIndex": " gs://YOUR_BUCKET/NORMAL.bai"
            },
            "tumorFinalBam": {
                "bam": " gs://YOUR_BUCKET/TUMOR.bam",
                "bamIndex": " gs://YOUR_BUCKET/TUMOR.bai"
            },
            "pairId": "TUMOR--NORMAL"
        }
    ],
    "normalSampleBamInfos": [
        {
            "sampleId": "NORMAL",
            "finalBam": {
                "bam": " gs://YOUR_BUCKET/NORMAL.bam",
                "bamIndex": " gs://YOUR_BUCKET/NORMAL.bai"
            }
        }
    ]
```

or they cam be taken from the `--pairs-file` optional columns `normalBam` and `tumorBam`

```
tumor,normal,pairId,tumorBam,normalBam
TUMOR,NORMAL,TUMOR--NORMAL,gs://YOUR_BUCKET/TUMOR.bam,gs://YOUR_BUCKET/NORMAL.bam
```

Use the `--dry-run` flag to validate the WDL, create an input JSON,
 validate the WDL with the input JSON and validate the uris to files
 
Use the `--skip-validate` flag to skid all validation of 
file uris. Disable with caution. Cromwell will launch
instances and run without checking. Test a small pairs
file to ensure all references exist and at least some
sample input files can be read by the current user.


### Run full pipeline from BAM
<a name="from_bam"></a>

Below are two examples of how to run the full pipeline starting from the BAM files as input.

```
# ===========================
# Run from BAMs with filled out json
# ===========================

program_dir="PATH_TO_REPO_PARENT/"
log_dir="PATH_TO/example/"
project_id="example_run_id"
cd ${log_dir}
bash ${program_dir}/wdl_port/run.sh \
--options ${program_dir}/wdl_port/example/options.json \
--wdl-file ${program_dir}/wdl_port/somatic_bam_wkf.wdl \
--url ${cromwell_server_url} \
--log-dir "${log_dir}" \
--project "${project_id}" \
--library WGS \
--genome Human_GRCh38_tcga \
--pairs-file ${program_dir}/wdl_port/example/tumorNormalPairs.csv \
--custom-inputs ${program_dir}/wdl_port/example/customInputs.json

# ===========================
# Run with BAM uris in csv template
# ===========================

program_dir="PATH_TO_REPO_PARENT/"
log_dir="PATH_TO/example/"
project_id="csv_run_id"
cd ${log_dir}
bash  ${program_dir}/wdl_port/run.sh \
--options  ${program_dir}/wdl_port/example/options.json \
--wdl-file  ${program_dir}/wdl_port/somatic_bam_wkf.wdl \
--url ${cromwell_server_url} \
--log-dir "${log_dir}" \
--project "${project_id}" \
--library WGS \
--genome Human_GRCh38_tcga \
--pairs-file  ${program_dir}/wdl_port/example/tumorNormalPairsWithBams.csv \
--custom-inputs  ${program_dir}/wdl_port/example/customInputsWithoutBams.json
```

### Run SNV, INDEL, SV and CNV calling from BAM
<a name="snv_indel_from_bam"></a>

Below is an example of how to run the SNV, INDEL, SV and CNV calling pipeline starting from the BAM files as input.

```
# ===========================
# Run SNV, INDEL, SV and CNV calling from BAM with uris in csv template
# ===========================
program_dir="PATH_TO_REPO_PARENT/"
log_dir="PATH_TO/example/"
project_id="csv_run_id"
cd ${log_dir}
bash ${program_dir}/wdl_port/run.sh \
--options ${program_dir}/wdl_port/example/options.json \
--wdl-file ${program_dir}/wdl_port/somatic_bam_gdc_wkf.wdl \
--url ${cromwell_server_url} \
--log-dir "${log_dir}" \
--project "${project_id}" \
--library WGS \
--genome Human_GRCh38_tcga \
--pairs-file ${program_dir}/wdl_port/example/tumorNormalPairsWithBams.csv \
--custom-inputs ${program_dir}/wdl_port/example/customInputsWithoutBams.json
```

### Run full pipeline from FASTQ
<a name="from_fastq"></a>

Below is an example of how to run the full pipeline starting from the FASTQ files as input.

```
# ===========================
# Run with FASTQ uris
# ===========================

program_dir="PATH_TO_REPO_PARENT/"
log_dir="PATH_TO/example/"
project_id="fastq_run_id"
cd ${log_dir}
bash ${program_dir}/wdl_port/run.sh \
--options ${program_dir}/wdl_port/example/options.json \
--wdl-file ${program_dir}/wdl_port/somatic_wkf.wdl \
--url ${cromwell_server_url} \
--log-dir "${log_dir}" \
--project "${project_id}" \
--library WGS \
--genome Human_GRCh38_tcga \
--pairs-file ${program_dir}/wdl_port/example/tumorNormalPairs.csv \
--custom-inputs ${program_dir}/wdl_port/example/customFastqInputs.json
```

