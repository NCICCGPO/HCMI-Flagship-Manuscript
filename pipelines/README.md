## HCMI Pipelines Index

This page lists the pipelines, tools, and commands used in the HCMI flagship analysis, organized by data modality. 

The `pipelines_launchers/` directory provides wrapper scripts for executing the corresponding pipelines, listed below. Each script follows the naming convention `run_pipeline_<pipeline_name>.<extension>`, matching the corresponding pipeline name.

---

### DNA Workflows (WGS/WES)

- **CN Concordance**: Scripts and workflows for assessing copy number (CN) concordance between models and their paired tumors.  
  Copy number concordance analysis *(Methods)*

- **Allele-Specific CNV + SV Inference (ReMixT)**: Inference of allele-specific CNVs and structural variants.
  https://github.com/amcpherson/remixt

- **ASCAT Segmentation**: Tumor–model CNV segmentation  
  The GDC pipeline for running ASCAT is described here: https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/DNA_Seq_WGS/

- **PURPLE Purity/Ploidy Calling**: Hartwig PURPLE pipeline  
  https://github.com/hartwigmedical/hmftools

- **Structural Variant Calls**  
  NYGC: https://github.com/nygenome/kancero
  
  Broad:
  
  WashU: https://github.com/ding-lab/SomaticSV
  
  MSKCC: https://github.com/amcpherson/remixt
  
  EMBL: https://github.com/hartwigmedical/hmftools

- **Copy Number Calls**  
  Broad: https://github.com/getzlab/ABSOLUTE
  
  WashU: https://github.com/mwyczalkowski/BICSEQ2.CWL
  
  MSKCC: https://github.com/amcpherson/remixt
  
  EMBL: https://github.com/hartwigmedical/hmftools

- **ecDNA Detection (AmpliconSuite)**: Inference of extrachromosomal DNA amplifications.  
  https://github.com/AmpliconSuite/AmpliconSuite-pipeline

- **Mutational Signature Analysis**  
  An implementation of the COSMIC signature framework: https://github.com/getzlab/SignatureAnalyzer

- **Driver Concordance Analysis**  
    1. [Setup Page](https://github.com/shahcompbio/driver-concordance/blob/main/README.md) - *Detailed instructions on required setup and running the pipeline.*
    2. Create Figures from Output Tables - *This can be done within the pipeline itself, or by running a downstream pipeline that generates figures from the intermediate outputs produced just prior to figure creation.*
        ```bash
        bash scripts/analyses/run_plotter_driver_concordance.sh
        ```
- **Compute purity/ploidy consensus**
  ```bash
  python scripts/compute_pp_consensus.py --input data/consensus_pp_for_github.txt --output_path '/path/for/output/'

- **SV consensus**: Computing consensus SV set using SV calls from individual centers. Modofied from PCAWG method. https://github.com/beroukhim-lab/hcmi_sv_consensus_public/tree/main
  ```bash
  #from linked repo, run on hcmi samples with
  bash scripts/run_sv_consensus_hcmi.bash
  ```
---

### RNA / Transcriptome Workflows

- **OncoMatch analysis**: Scripts to perform protein activity-based OncoMatch analysis on the HCMI collection.

- **Celligner analysis**: Pipeline to align HCMI RNA expression profiles with TCGA/TARGET and CCLE datasets using the Celligner framework as described in Warren et al., *Nature Communications* 2021 (https://doi.org/10.1038/s41467-020-20294-x) and https://github.com/broadinstitute/celligner.

- **MultiClass Pairs Based Distance Analysis**: *Tumor matched model distance calculation based on multiclass pair prediction scores*

- **Latent Transcription Factor Distance Analysis**: *Sample biological distance calculation based on transcription factors.*
    1. [Setup Page](https://github.com/jordan2lee/latent-tf-tumor-growth-distances/blob/main/doc/requirements.md) - *Detailed instructions on required set up.*
    2. Run Pipeline - *Runs all cancers for an inter-cohort analysis*
        ```bash
        bash scripts/analyses/run_pipeline_latenttf_dist.sh
        ```
        Additional information on [Analysis Page](https://github.com/jordan2lee/latent-tf-tumor-growth-distances).

- **Euclidean distance analysis**: Sample biological distance calculation and analysis. 
    1. [Setup Page](https://github.com/jordan2lee/euclidean-tumor-growth-distance/blob/main/doc/requirements.md) - Detailed instructions on required set up.
    2. Run Pipeline - Specify a single cancer cohort or use `ALL` to run all cancer cohorts (all options found [here](https://github.com/jordan2lee/euclidean-tumor-growth-distance/blob/main/doc/cohort_options.md)). Second specify name of folder to house intermediate files in pipeline submodule.
        ```bash
        bash scripts/analyses/run_pipeline_euc_dist.sh PAAD data/prep
        ```
        Additional information on [Analysis Page](https://github.com/jordan2lee/euclidean-tumor-growth-distance).

---

### Epigenetics / DNA Methylation Workflows

- **TMP Toolkit Subtype Classification**: Labeling samples with TMP subtypes.
    + [Setup Page](https://github.com/jordan2lee/classify-lab-models-and-tumors/blob/main/doc/requirements.md) - *Detailed instructions on required set up.*
    + Run Pre-processing Pipeline - *Specify a single cancer cohort or use ALL to run all cancer cohorts (all options found [here](https://github.com/jordan2lee/classify-lab-models-and-tumors/blob/main/doc/cohort_options.md)). Second specify name of folder to house intermediate files in pipeline submodule.*
        ```bash
        bash scripts/analyses/run_pipeline_preprocess_tmp_subtype.sh PAAD data/prep
        ```
        Additional information on [Analysis Page](https://github.com/jordan2lee/classify-lab-models-and-tumors/).
    
    + Run Classification Pipeline - *Specify a single cancer cohort*
        ```bash
        # Run Gene Expression based Classifier
        bash run_pipeline_get_tmp_subtype.sh PAAD data/prep GEXP
        ```
        or
        ```bash
        # Run DNA Methylation based Classifier
        bash run_pipeline_get_tmp_subtype.sh PAAD data/prep METH
        ```
