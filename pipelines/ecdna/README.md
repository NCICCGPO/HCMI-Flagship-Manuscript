# GDAN HCMI Project: ecDNA detection pipeline

## Overview
This directory contains the code necessary for generating ecDNA calls for the HCMI cohort, as well as the two ecDNA-related Extended Data figures in the associated manuscript. This pipeline calls ecDNAs and other focal amplifications via [AmpliconSuite-pipeline](https://github.com/AmpliconSuite/AmpliconSuite-pipeline), using the GroupedAnalysis module for tumor/model pairs where possible. ecDNAs are then compared between tumors and models to generate a set of concordant calls, and ecDNAs that are called in both samples but only marked PASS in one are "rescued" back into the callset. 


## Dependencies
* A Slurm HPC environment
* [AmpliconSuite-pipeline](https://github.com/AmpliconSuite/AmpliconSuite-pipeline)
* R 4.0.0 or greater
* The following R packages
	* optparse
    * reshape2
	* stringr
	* gUtils
	* ggplot2
	* ggbeeswarm


## Running the pipeline  
Running the pipeline requires that you obtain the JaBbA copy number calls and fragCounter coverage profiles, as well as the BAMs and indices for the HCMI cohort. After you have these, you can edit `PROJECT_DIR` to point to the location where the output directory will be made. You will also need to create a BAM mapping file (see example `data/nygc_bam_map.csv`) to map between GDC IDs and BAM locations. After the data is in place, simply run:

```{bash}
bash scripts/run.sh
```