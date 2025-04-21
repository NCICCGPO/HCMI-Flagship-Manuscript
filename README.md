# HCMI-Flagship-Manuscript  

This repository contains software and analysis code used for the publication:  
**A Compendium of Cancer Organoid Models for Diverse Cancer Types**  

The repository is organized into multiple Git submodules, each pointing to individual repositories and the version of the code used for the analysis.  

## Contents  

This repository currently includes the following components:  

- **CNV Concordance Analysis**: Scripts and workflows for assessing copy number variation (CNV) concordance between models and their paired tumors.  
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

## Cloning This Repository  

To properly check out this repository and ensure all submodules are included, use the following command:  

```bash
git clone --recurse-submodules https://github.com/NCICCGPO/HCMI-Flagship-Manuscript.git
```  

This will automatically retrieve all submodules required for the analysis.  
