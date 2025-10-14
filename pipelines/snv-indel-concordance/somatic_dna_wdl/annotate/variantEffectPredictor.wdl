version 1.0

import "../wdl_structs.wdl"

task vepPublicSvnIndel {
    input {
        IndexedReference vepFastaReference
        String pairName
        IndexedVcf unannotatedVcf

        # Somatic
        IndexedVcf cosmicCoding
        IndexedVcf cosmicNoncoding
        IndexedVcf gnomadGenomes
        IndexedVcf gnomadExomes

        # Public
        File vepCache
        File annotations
        File plugins
        File dbNSFPReplacementLogic
        File MaxEntScan
        String vepGenomeBuild
        IndexedTable dbNSFP4
        IndexedTable dbscSNV

        # Public
        IndexedVcf deepIntronicsVcf
        IndexedVcf clinvarIntronicsVcf

        String vcfAnnotatedVepPath = "~{pairName}.v7.vep.annotated.vcf"
        Int threads = 4
        Int memoryGb = 7
        Int diskSize

    }

    command {
        set -e -o pipefail

        # NOTE:task will not work with any other genome build as is
        # because of this section

        mkdir -p ensembl_vep
        cd ensembl_vep
        tar -xzvf ~{vepCache}
        tar -xzvf ~{annotations}
        tar -xzvf ~{plugins}
        tar -xzvf ~{MaxEntScan}
        
        cd ../

        /opt/vep/src/ensembl-vep/vep \
        --fork ~{threads} \
        --buffer_size 50000 \
        --format vcf \
        --no_stats \
        --no_escape \
        --offline \
        --assembly ~{vepGenomeBuild} \
        --cache \
        --dir_cache ensembl_vep \
        --refseq \
        --max_af \
        --af \
        --af_1kg \
        --af_gnomad \
        --exclude_predicted \
        --fasta ~{vepFastaReference.fasta} \
        --symbol \
        --hgvs \
        --check_existing \
        --vcf \
        --pick_allele_gene \
        --dir_plugins ensembl_vep/Plugins \
        --plugin dbscSNV,~{dbscSNV.table} \
        --plugin MaxEntScan,ensembl_vep/fordownload \
        --plugin dbNSFP,~{dbNSFP4.table},~{dbNSFPReplacementLogic},REVEL_score,SIFT_pred,SIFT4G_pred,LRT_pred,MutationTaster_pred,MutationAssessor_pred,FATHMM_pred,PROVEAN_pred,MetaSVM_pred,PrimateAI_pred,fathmm-MKL_coding_pred,GERP++_RS,phyloP100way_vertebrate,CADD_phred,Polyphen2_HVAR_pred \
        --custom ensembl_vep/annotations/clinvar.vep.vcf.gz,CLN_Overlap,vcf,overlap,0,CLIN_ID,CLNSIG,CLNREVSTAT,CLNDN \
        --custom ensembl_vep/annotations/clinvar.vep.vcf.gz,CLN_Exact,vcf,exact,0,CLIN_ID,CLNSIG,CLNREVSTAT,CLNDN \
        --custom ~{gnomadGenomes.vcf},GnomadGenomes,vcf,exact,0,AF_fin_XX,AF_nfe_XX,AF_oth_XX,AF_sas_XX,AF_ami_XX,AF_asj_XX,AF_mid_XX,AF_eas_XX,AF_amr_XX,AF_afr_XX,AF_XX,AF_fin_XY,AF_nfe_XY,AF_oth_XY,AF_sas_XY,AF_ami_XY,AF_asj_XY,AF_mid_XY,AF_eas_XY,AF_amr_XY,AF_afr_XY,AF_XY,AF_fin,AF_nfe,AF_oth,AF_sas,AF_ami,AF_asj,AF_mid,AF_eas,AF_amr,AF_afr,AF,AF_non_cancer_fin_XX,AF_non_cancer_nfe_XX,AF_non_cancer_oth_XX,AF_non_cancer_sas_XX,AF_non_cancer_ami_XX,AF_non_cancer_asj_XX,AF_non_cancer_mid_XX,AF_non_cancer_eas_XX,AF_non_cancer_amr_XX,AF_non_cancer_afr_XX,AF_non_cancer_XX,AF_non_cancer_fin_XY,AF_non_cancer_nfe_XY,AF_non_cancer_oth_XY,AF_non_cancer_sas_XY,AF_non_cancer_ami_XY,AF_non_cancer_asj_XY,AF_non_cancer_mid_XY,AF_non_cancer_eas_XY,AF_non_cancer_amr_XY,AF_non_cancer_afr_XY,AF_non_cancer_XY,AF_non_cancer_fin,AF_non_cancer_nfe,AF_non_cancer_oth,AF_non_cancer_sas,AF_non_cancer_ami,AF_non_cancer_asj,AF_non_cancer_mid,AF_non_cancer_eas,AF_non_cancer_amr,AF_non_cancer_afr,AF_non_cancer,AN,AN_non_cancer,nhomalt,nhomalt_non_cancer \
        --custom ~{gnomadExomes.vcf},GnomadExomes,vcf,exact,0,AF_fin_female,AF_nfe_female,AF_oth_female,AF_sas_female,AF_asj_female,AF_eas_female,AF_amr_female,AF_afr_female,AF_female,AF_fin_male,AF_nfe_male,AF_oth_male,AF_sas_male,AF_asj_male,AF_eas_male,AF_amr_male,AF_afr_male,AF_male,AF_fin,AF_nfe,AF_oth,AF_sas,AF_asj,AF_eas,AF_amr,AF_afr,AF,non_cancer_AF_fin_female,non_cancer_AF_nfe_female,non_cancer_AF_oth_female,non_cancer_AF_sas_female,non_cancer_AF_asj_female,non_cancer_AF_eas_female,non_cancer_AF_amr_female,non_cancer_AF_afr_female,non_cancer_AF_female,non_cancer_AF_fin_male,non_cancer_AF_nfe_male,non_cancer_AF_oth_male,non_cancer_AF_sas_male,non_cancer_AF_asj_male,non_cancer_AF_eas_male,non_cancer_AF_amr_male,non_cancer_AF_afr_male,non_cancer_AF_male,non_cancer_AF_fin,non_cancer_AF_nfe,non_cancer_AF_oth,non_cancer_AF_sas,non_cancer_AF_asj,non_cancer_AF_eas,non_cancer_AF_amr,non_cancer_AF_afr,non_cancer_AF,AN,non_cancer_AN,nhomalt,non_cancer_nhomalt \
        --custom ~{cosmicCoding.vcf},CosmicCoding,vcf,exact,0,GENE,TRANSCRIPT,LEGACY_ID,SAMPLE_COUNT,TIER,IS_CANONICAL,HGVSG,AA,CDS,HGVSP \
        --custom ~{cosmicNoncoding.vcf},CosmicNonCoding,vcf,exact,0,GENE,TRANSCRIPT,LEGACY_ID,SAMPLE_COUNT,TIER,IS_CANONICAL,HGVSG \
        --custom ~{deepIntronicsVcf.vcf},INTRONIC,vcf,exact,0,INTRONIC \
        --custom ~{clinvarIntronicsVcf.vcf},CLINVAR_INTRONIC,vcf,exact,0,INTRONIC \
        --custom ensembl_vep/annotations/spliceai_scores.hg38.sorted.vcf.gz,SPLICEAI,vcf,exact,0,DS_AG,DS_AL,DS_DG,DS_DL \
        --custom ensembl_vep/annotations/pli_hg38.vcf.gz,PLI,vcf,overlap,0,pLI,mis_z \
        --custom ensembl_vep/annotations/dials_genes_b38.vcf.gz,DIALS,vcf,overlap,0,DIALS_GENE \
        --input_file ~{unannotatedVcf.vcf} \
        --output_file ~{vcfAnnotatedVepPath}

    }

    runtime {
        mem: memoryGb + "G"
        cpus: threads
        cpu : threads
        memory : memoryGb + "GB"
        disks: "local-disk " + diskSize + " HDD"
        docker: "gcr.io/nygc-public/ensemblorg/ensembl-vep@sha256:7c99fadb332d65fc913844ce93aaee89dac5516afb68179304f1973a9e270196"
        runtime_minutes: "500"
    }

    output {
        File vcfAnnotatedVep = "~{vcfAnnotatedVepPath}"
    }
}
