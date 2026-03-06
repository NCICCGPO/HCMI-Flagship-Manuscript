- **CN Concordance**: Scripts and workflows for assessing copy number (CN) concordance between models and their paired tumors.  
  Copy number concordance analysis (see methods for detailed description on how ABSOLUTE forcecalling was performed).
  A basic procedure for ABSOLUTE forcecalling:
  ```
  bash ~pipelines/consensus_pipelines/compute_cnv_consensus.bash
  ```

- **Driver Concordance Analysis**  
    1. [Setup Page](https://github.com/shahcompbio/driver-concordance/blob/main/README.md) - *Detailed instructions on required setup and running the pipeline.*
    2. Create Figures from Output Tables - *This can be done within the pipeline itself, or by running a downstream pipeline that generates figures from the intermediate outputs produced just prior to figure creation.*
        ```bash
        bash scripts/analyses/run_plotter_driver_concordance.sh
        ```
- **Compute purity/ploidy consensus**
  ```bash
  python ~pipelines/consensus_pipelines/compute_pp_consensus.py --input data/consensus_pp_for_github.txt --output_path '/path/for/output/'

- **SV consensus**: Computing consensus SV set using SV calls from individual centers. Modofied from PCAWG method. https://github.com/beroukhim-lab/hcmi_sv_consensus_public/tree/main
  ```bash
  #from linked repo, run on hcmi samples with
  bash scripts/run_sv_consensus_hcmi.bash
  ```
