# HCMI-Flagship-Manuscript  

This repository contains software and analysis code used for the publication:  
**A Compendium of Cancer Organoid Models for Diverse Cancer Types**  

The repository is organized into multiple Git submodules, each pointing to individual repositories and the version of the code used for the analysis.  

## Contents  

This repository currently includes the following components:  

- **CNV Concordance Analysis**: Scripts and workflows for assessing copy number variation (CNV) concordance between models and their paired tumors.  
    1. [Setup Page](https://github.com/shahcompbio/driver-concordance/blob/main/README.md) - *Detailed instructions on required setup and running the pipeline.*
    2. Create Figures from Output Tables - *This can be done within the pipeline itself, or by running a downstream pipeline that generates figures from the intermediate outputs produced just prior to figure creation.*
        ```bash
        bash scripts/analyses/run_plotter_driver_concordance.sh
        ```

- **OncoMatch analysis**: Scripts to perform protein activity-based OncoMatch analysis on the HCMI collection.
- **Single nuclei analysis**: Demo codes highlighting snRNA-seq data analysis in selected tumor-model pairs.
- **HCMI Explorer Suite**: A shiny app to explore treatment timelines of the HCMI collection.  
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


## Cloning This Repository  

To properly check out this repository and ensure all submodules are included, use the following command:  

```bash
git clone --recurse-submodules https://github.com/NCICCGPO/HCMI-Flagship-Manuscript.git
```  

This will automatically retrieve all submodules required for the analysis.  
