# HCMI-Flagship-Manuscript  

This repository contains software and analysis code used for the publication:  
**A Compendium of Cancer Organoid Models for Diverse Cancer Types**  

The repository is organized into multiple Git submodules, each pointing to individual repositories and the version of the code used for the analysis.  

## Contents  

This repository currently includes the following components:  

- **CNV Concordance Analysis** – Scripts and workflows for assessing copy number variation (CNV) concordance between models and their paired tumors.  
- **HCMI Explorer Suite** – A shiny app to explore treatment timelines of the HCMI collection.  
- **OncoMatch analysis** - Scripts to perform protein activity-based OncoMatch analysis on the HCMI collection.
- **Single nuclei analysis** - Demo codes highlighting snRNA-seq data analysis in selected tumor-model pairs.
  
## Cloning This Repository  

To properly check out this repository and ensure all submodules are included, use the following command:  

```bash
git clone --recurse-submodules https://github.com/NCICCGPO/HCMI-Flagship-Manuscript.git
```  

This will automatically retrieve all submodules required for the analysis.  
