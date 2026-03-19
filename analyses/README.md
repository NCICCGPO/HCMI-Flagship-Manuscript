## HCMI Analyses Index

This page provides an overview of the analysis workflows, notebooks, and figure-generation codes used in the HCMI flagship analysis.

---

### notebooks

This directory contains analytical workflow and notebooks (.ipynb, .html, .qmd file extensions) used for data processing, downstream analysis and figure generation. 

**Fig.3a**: `Epigenetic_Concordance_Fig.3a.html`

This script calculates the Pearson correlation between models and parent tumors based on their cancer-associated CpG hypermethylation profiles; these results are used to generate Figure 3a.

**Fig4ab**: `Fig4ab.ipynb`

This notebook generates Figure 4a–b by summarizing the distribution of disease types across HCMI models and CCLE cell lines. It aggregates RNA-seq metadata by disease and model type, calculates relative representation, and produces the comparison plots shown in Fig. 4a and the summary of unique HCMI-only disease subtypes shown in Fig. 4b.

**Fig.4c**: `Fig.4c_NMF.html` + `Fig.4c_UMAP.html`

These scripts perform dimension reduction via non-negative matrix factorization (NMF) on combined HCMI, TCGA, and TARGET cohorts using CpG sites with cancer-associated DNA hypermethylation. A UMAP is then generated based on the NMF factors for tumor types represented across both HCMI and the TCGA/TARGET projects

**Fig4d**: `Fig4d.ipynb`

This notebook generates Figure 4d using UMAP coordinates from Celligner outputs generated with the Celligner pipeline. It visualizes TCGA/TARGET tumors, HCMI tumors, and HCMI models in a shared transcriptomic embedding space, with colors denoting cancer type and symbols distinguishing sample types.

**Fig.4e**: `Fig.4e.html`

This script generates a heatmap of cancer-associated DNA hypermethylation profiles for HCMI and TCGA samples from four major cancer cohorts.

**Fig.4f**: `HCMI-Flagship-Manuscript/analyses/notebooks/Fig4f_heatmap.html` 

This script visualizes the heatmap of protein activities for top 4 cancer cohorts (COAD/READ, STAD/ESCA/, GBM/LGG, PAAD)

**Fig.4g** + **Fig.E5f**: `HCMI-Flagship-Manuscript/analyses/notebooks/Fig4g-ExtendedFig5f_violinplots-HCMI-vs-CCLE.html`

This script visualizes the distribution of OncoMatch NES scores (representing TCGA tumor–HCMI/CCLE model similarity) using violin plots for common cancer types in the HCMI and CCLE datasets.

**Fig.E5c**: `Extended_Data_Fig.5c.html`

This script generates a UMAP using the same NMF factors as Figure 4c, covering all tumor types represented across the HCMI, TCGA, and TARGET projects.

**Fig. 5** + **Figs. E6, E7, E8**: `HCMI-single-nuclei` submodule

`HCMI-single-nuclei` provides a standalone, conventiently organized, self-contained submodule that implements the complete snRNA-seq analytical workflow. It includes all notebooks required for preprocessing, downstream analysis, and figure generation of single-nuclei RNA-seq data and a dedicated README.

**Extended_data_5d**: `Extended_data_5d.ipynb`

This notebook generates Extended Data Figure 5d using  UMAP coordinates from the Celligner pipeline to further assess transcriptomic relationships between tumors and models. This visualization includes the CCLE Cell Lines in addition to the TCGA/TARGET tumors, HCMI tumors, and HCMI models shown in Figure 4d.



---

### scripts

This directory contains executable scripts supporting figure generation for Figure 2 and Extended Data Figure 5. 

**Fig.3e** + **Fig.E4d**: `HCMI-Flagship-Manuscript/analyses/scripts/Fig3e-ExtendedFig4d-blob_generator.R` 

Visualization of the top ten candidate MR proteins and expression of target genes for a representative high-fidelity GBM model (HCM-BROD-0689-C71). Extended Fig. 4d for an example of a mismatch.
