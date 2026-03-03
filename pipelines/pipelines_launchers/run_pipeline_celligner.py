"""
HCMI Celligner Alignment Pipeline
-------------------------------

This script performs transcriptomic alignment of the HCMI tumor and model datasets with TCGA/TARGET and CCLE
using the Celligner algorithm, following the protocol described in:

Warren et al., *Nature Communications* 2021 (https://doi.org/10.1038/s41467-020-20294-x) and https://github.com/broadinstitute/celligner

All input matrices are log2(TPM + 1)-transformed gene expression profiles, filtered to retain
only protein-coding genes shared across datasets. 

"""

import celligner
import pandas as pd
import numpy as np


# ----------------------------------------
# Load Input Datasets (log2(TPM+1) format)
# ----------------------------------------

#TCGA_expression_input--> TARGET program samples (n = 784) and TCGA expression data (n = 9,806) were obtained from the Xena browser (https://xenabrowser.net)
#CCLE_expression_input--> CCLE dataset for 1377 samples were downloaded from the DepMap portal Public 19Q4 release
#HCMI_expression_input--> HCMI RNA expression data is available on GDC publication page


# Create Celligner object and fit + transform the initial reference (CCLE) and target (TCGA) expression datasets
mycelligner = celligner.Celligner()
mycelligner.fit(CCLE_expression_input)
mycelligner.transform(TCGA_expression_input)

# Compute UMAP, clusters and tumor - model distance
mycelligner.computeMetricsForOutput()

# Make new reference from CCLE and TCGA and then bring in the HCMI dataset
mycelligner.makeNewReference()
mycelligner.ref_clusters=mycelligner.output_clusters
mycelligner.mnn_kwargs.update({"k1":20, "k2":50}) 
mycelligner.cpca_ncomp=3
mycelligner.transform(HCMI_expression_input, compute_cPCs=True)

# Compute UMAP, and Celligner pairwise euclidean distances at high dimensional space
mycelligner.computeMetricsForOutput()

# ----------------------------------------
# Outputs
# ----------------------------------------

# All Celligner outputs used in the manuscript are deposited in the GDC publication page:
# - Supplementary Table 11: RNA UMAP coordinates sheet of all samples
# - Supplementary Table 12: Pairwise Celligner distances (70-PC Euclidean space) of HCMI models and tumors with lineage assignments
# - Celligner-aligned expression matrix is provided in GDC publication page: Supplementary_Table_Celligner_aligned_data.gz
# - Celligner loadings are provided in GDC publication page: Supplementary_Table_Celligner_cPC_loadings.xlsx
# - Celligner pairwise Euclidean distances in 70 PC space for all samples are provided in GDC publication page: Supplementary_Table_Celligner_pairwise_euclidean_distances_in_70PC.csv

# Reference:
# Warren A, et al. Celligner: A computational method to align tumor and cell line transcriptomes. Nat Commun. 2021.
# Main Celligner GitHub repository for Celligner source code and functions: https://github.com/broadinstitute/celligner

