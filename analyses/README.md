## HCMI Analyses Index

This page provides an overview of the analysis workflows, notebooks, and figure-generation codes used in the HCMI flagship analysis.

---

### notebooks 
This directory contains analytical workflow and notebooks (.ipynb, .html, .qmd file extensions) used for data processing, downstream analysis and figure generation. 


-**Fig.4f**: `HCMI-Flagship-Manuscript/analyses/notebooks/Fig4f_heatmap.html` 

This script visualizes the heatmap of protein activities for top 4 cancer cohorts (COAD/READ, STAD/ESCA/, GBM/LGG, PAAD)


-**Fig.4g** + **Fig.E5f**: `HCMI-Flagship-Manuscript/analyses/notebooks/Fig4g-ExtendedFig5f_violinplots-HCMI-vs-CCLE.html`

This script visualizes the distribution of OncoMatch NES scores (representing TCGA tumor–HCMI/CCLE model similarity) using violin plots for common cancer types in the HCMI and CCLE datasets.

-**Fig. 5** + **Figs. E6, E7, E8**: `HCMI-single-nuclei` submodule

`HCMI-single-nuclei` provides a standalone, conventiently organized, self-contained submodule that implements the complete snRNA-seq analytical workflow. It includes all notebooks required for preprocessing, downstream analysis, and figure generation of single-nuclei RNA-seq data and a dedicated README.

---

-### scripts
This directory contains executable scripts supporting figure generation for Figure 2 and Extended Data Figure 5. 

-**Fig.3e** + **Fig.E4d**: `HCMI-Flagship-Manuscript/analyses/scripts/Fig3e-ExtendedFig4d-blob_generator.R` 

Visualization of the top ten candidate MR proteins and expression of target genes for a representative high-fidelity GBM model (HCM-BROD-0689-C71). Extended Fig. 4d for an example of a mismatch.



