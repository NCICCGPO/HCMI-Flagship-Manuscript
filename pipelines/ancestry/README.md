# NYGC Genetic ancestry pipeline 

Ancestry proportion was determined by ADMIXTURE (v1.3.0), which used a maximum likelihood-based method to estimate the proportion of reference-population ancestries in a sample. To do this we genotyped reference markers that we generated from 1,964 unrelated 1000 Genomes project samples directly on the whole genome samples using GATK pileup (v3.4.0). We excluded individuals from populations MXL (Mexican Ancestry from Los Angeles, USA), ACB (African Caribbean in Barbados), and ASW (African Ancestry in the Southwest US) from the reference due to their being putatively admixed. We further filtered the reference by using only SNP markers with a minimum minor allele frequency (MAF) of 0.01 overall and 0.05 in at least one 1000Genomes continental population. Variants were additionally linkage disequilibrium (LD)-pruned using PLINK (v1.9) with a window size of 500kb, a step size of 250kb and an r2 threshold of 0.2. The analysis resulted in a proportional breakdown of each sample into five continental populations (AFR, AMR, EAS, EUR, SAS) and 23 populations. We then categorized patients by the continental population of highest proportion. 


**VERY IMPORTANT: ADMIXTURE fails at some step if there is an underscore "_" in the sample name. If there is an underscore in the sample name, change that to some other character in the -sample field.**

**Required resources**
Please download the following resources and modify the json config file before running
```
GRCh38_full_analysis_set_plus_decoy_hla.fa: gs://nygc-resources-public/GRCh38_full_analysis_set_plus_decoy_hla/internal/GRCh38_full_analysis_set_plus_decoy_hla.fa
dbSNP b147_b38 00-All.vcf.gz: gs://nygc-resources-public/dbSNP/b147_b38/00-All.vcf.gz 
```
````
resources/cm/markers/path_to_reference_and_markers_for_genotyping.json
```
**Example usage :**

```
sbatch scripts/ancestry_admixture.sh \
 -bam [SAMPLE_BAM] \
 -sample [SAMPLE_NAME] \
 -outdir [OUTDIR] \
 -seq_type wgs \
 -subpops
 ```
