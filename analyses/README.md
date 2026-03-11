## HCMI Analyses Index

This page provides an overview of the analysis workflows, notebooks, and figure-generation codes used in the HCMI flagship analysis.

---

### notebooks 
This directory contains analytical workflow and notebooks (.ipynb, .html, .qmd file extensions) used for data processing, downstream analysis and figure generation. `HCMI-single-nuclei` provides a standalone, conventiently organized, self-contained submodule that implements the complete snRNA-seq analytical workflow. It includes all notebooks required for preprocessing, downstream analysis, and figure generation of single-nuclei RNA-seq data.

**Fig.4f**: `HCMI-Flagship-Manuscript/analyses/notebooks/Fig4f_heatmap.html` 
This script visualizes the heatmap of protein activities for top 4 cancer cohorts (COAD/READ, STAD/ESCA/, GBM/LGG, PAAD)


**Fig.4g** + **Fig.E5f**: `HCMI-Flagship-Manuscript/analyses/notebooks/Fig4g-ExtendedFig5f_violinplots-HCMI-vs-CCLE.html`
This script visualizes the distribution of OncoMatch NES scores (representing TCGA tumor–HCMI/CCLE model similarity) using violin plots for common cancer types in the HCMI and CCLE datasets.

---

### scripts
This directory contains executable scripts supporting figure generation for Figure 2 and Extended Data Figure 5. 



