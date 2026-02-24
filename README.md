# HCMI-Flagship-Manuscript  

This repository contains software and analysis code used for the publication:  
**A Compendium of Next-Generation Patient-Derived Models for Diverse Cancer Types**  

### About HCMI

The Human Cancer Models Initiative (HCMI) is an international consortium generating next-generation patient-derived cancer models with matched molecular and clinical data. The initiative expands disease representation and provides standardized, deeply characterized tumor-derived models to support research and translational discovery.

More [about HCMI](https://www.cancer.gov/ccg/research/functional-genomics/hcmi)

---

### HCMI Data Access

Public and controlled-access data, including WGS/WES, RNA-seq, DNA methylation, WSIs, and clinical metadata, are available from the [NCI Genomic Data Commons (GDC)](https://portal.gdc.cancer.gov/projects/HCMI-CMDC)

---

### HCMI Analysis Pipelines & Code

This repository includes analysis scripts and pipeline references used to generate the results presented in the HCMI flagship manuscript.

- **Analysis Pipelines** — 
Computational workflows for tumor–model concordance and biological distance metrics.  
Pipelines are included as submodules and documented at: `pipelines/`

- **Analysis Scripts** — 
Data preprocessing, analytical workflows, and supporting scripts used for data processing and figure generation: `scripts/analyses/`. This directory also contains the `HCMI-single-nuclei` submodule, which includes all Jupyter notebooks for snRNA-seq analyses.

- **Visualization & Manuscript Figure Scripts** — 
Scripts used to generate figures included in the manuscript: `scripts/figures/`

---

### HCMI Platforms & Interactive Tools

The HCMI ecosystem includes multiple platforms for exploring, accessing, and requesting patient-derived models:

- **[HCMI Searchable Catalogue](https://hcmi-searchable-catalog.nci.nih.gov/)** — primary portal to browse available HCMI models and associated metadata  

- **HCMI Explorer Suite** — interactive viewer for clinical timelines, and genomic annotations. Platform repository submoduled under: `platforms/`

- **[cBioPortal (HCMI Collection)](https://www.cbioportal.org/)** — integrated genomic visualization for HCMI models  

- **[HCMI page @ ATCC](https://www.atcc.org/hcmi)** — models, culturing protocols, and ordering information available through ATCC  

---

### Cloning This Repository  

To properly check out this repository and ensure all submodules are included, use the following command:  

```bash
git clone --recurse-submodules https://github.com/NCICCGPO/HCMI-Flagship-Manuscript.git
```  

This will automatically retrieve all submodules required for the analysis.  

---


