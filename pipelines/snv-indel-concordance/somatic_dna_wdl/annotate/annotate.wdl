version 1.0

import "../wdl_structs.wdl"

# General tasks

task MainCenterVcf {
    input {
        Int memoryGb = 2
        String pairName
        String mainVcfPath = "~{pairName}.snv.indel.union.vcf"
        String mainIntersectionVcfPath = "~{pairName}.snv.indel.intersection.vcf"
        String mainTmpVcfPath = "~{pairName}.snv.indel.union.tmp.vcf"
        File vcfAnnotated
        Int diskSize = ceil( size(vcfAnnotated, "GB") * 2) + 5
        # File multiCenterAnnotation = "pipelines/snv-indel-concordance/somatic_dna_tools/multi_center_annotation.py"
        File multiCenterAnnotation = "gs://nygc-comp-s-fd4e-input/somatic/hcmi/tools/multi_center_annotation.py"
    }

    command {
        bcftools annotate \
        --remove "INFO/nygc_CSQ,INFO/washu_CSQ,INFO/broad_CSQ,INFO/washuExome_CSQ,INFO/broadExome_CSQ,INFO/washuWgs_CSQ,INFO/broadWgs_CSQ" \
        ~{vcfAnnotated} \
        > ~{mainTmpVcfPath}
        
        python \
        ~{multiCenterAnnotation} \
        --vcf ~{mainTmpVcfPath} \
        --multi-caller-out ~{mainIntersectionVcfPath} \
        --out ~{mainVcfPath}
    }

    output {
        File mainVcf = "~{mainVcfPath}"
        File mainIntersectionVcf = "~{mainIntersectionVcfPath}"
    }

    runtime {
        mem: memoryGb + "G"
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-comp-s-fd4e/gcs_htslib_suite@sha256:47eba58683641905a58f31c05605c953dc1d888466e8d434a6ceff47b76df03e"
        runtime_minutes: "1440"
    }
}

task MainExonicCenterVcf {
    input {
        Int memoryGb = 2
        String pairName
        String mainVcfPath = "~{pairName}.snv.indel.union.vcf"
        String mainIntersectionVcfPath = "~{pairName}.snv.indel.intersection.vcf"
        String mainTmpVcfPath = "~{pairName}.snv.indel.union.tmp.vcf"
        File vcfAnnotated
        IndexedTable exonicRegionBedAnnotationsGz
        String tmpVcfHeaderPath = "vcf_header_lines.hdr"
        Int diskSize = ceil( size(vcfAnnotated, "GB") * 2) + 5
        #File multiCenterAnnotation = "pipelines/snv-indel-concordance/somatic_dna_tools/multi_center_annotation.py"
        File multiCenterAnnotation = "gs://nygc-comp-s-fd4e-input/somatic/hcmi/tools/multi_center_annotation.py"
    }

    command {
        set -e -o pipefail
        echo \
        '##INFO=<ID=ExonicRegion,Number=0,Type=Flag,Description="Variant in a Exonic region">' \
        > ~{tmpVcfHeaderPath}

        bcftools annotate \
        --remove "INFO/nygc_CSQ,INFO/washu_CSQ,INFO/broad_CSQ,INFO/washuExome_CSQ,INFO/broadExome_CSQ,INFO/washuWgs_CSQ,INFO/broadWgs_CSQ" \
        --header-lines ~{tmpVcfHeaderPath} \
        --columns CHROM,FROM,TO,ExonicRegion \
        --annotations ~{exonicRegionBedAnnotationsGz.table} \
        --min-overlap :1 \
        ~{vcfAnnotated} \
        > ~{mainTmpVcfPath}
        
        python \
        ~{multiCenterAnnotation} \
        --vcf ~{mainTmpVcfPath} \
        --multi-caller-out ~{mainIntersectionVcfPath} \
        --out ~{mainVcfPath} \
        --exome-wgs
    }

    output {
        File mainVcf = "~{mainVcfPath}"
        File mainIntersectionVcf = "~{mainIntersectionVcfPath}"
    }

    runtime {
        mem: memoryGb + "G"
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-comp-s-fd4e/gcs_htslib_suite@sha256:47eba58683641905a58f31c05605c953dc1d888466e8d434a6ceff47b76df03e"
        runtime_minutes: "1440"
    }
}



task GermFilterVcf {
    input {
        String pairName
        String chrom
        String filteredOutFilePath = "~{pairName}.final.v7.filtered.~{chrom}.vcf"
        File ponOutFile
        Int memoryGb = 60
        Int diskSize = 16
    }

    command {
        python \
        /germ_filter_vcf.py \
        ~{ponOutFile} \
        ~{filteredOutFilePath}
    }

    output {
        File filteredOutFile = "~{filteredOutFilePath}"
    }

    runtime {
        mem: memoryGb + "G"
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-comp-s-fd4e/gcs_htslib_suite@sha256:47eba58683641905a58f31c05605c953dc1d888466e8d434a6ceff47b76df03e"
        runtime_minutes: "500"
    }
}

task CustomAnnotateCenter {
    input {
        String pairName
        File vcfAnnotatedVep
        String vcfCsqRenamedPath = "~{pairName}.snv.indel.renamed.annotated.vcf"
        String vcfSupplementalPath = "~{pairName}.snv.indel.centers.annotated.vcf"
        Int memoryGb = 2
        Int diskSize = ceil( size(vcfAnnotatedVep, "GB") * 2) + 5
        # File germFilterCenterVcf = "pipelines/snv-indel-concordance/somatic_dna_tools/germ_filter_center_vcf.py"
        File germFilterCenterVcf = "gs://nygc-comp-s-fd4e-input/somatic/hcmi/tools/germ_filter_center_vcf.py"

    }

    command {
    set -e -o pipefail
    # RenameCsqVcf
        python3 \
        /rename_csq_vcf.py \
        ~{vcfAnnotatedVep} \
        ~{vcfCsqRenamedPath}
    # GermFilterVcf
        python3 \
        ~{germFilterCenterVcf} \
        ~{vcfCsqRenamedPath} \
        ~{vcfSupplementalPath}
    }
    
    output {
        File vcfCsqRenamed = "~{vcfSupplementalPath}"
    }
    
    runtime {
        mem: memoryGb + "G"
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-comp-s-fd4e/gcs_htslib_suite@sha256:47eba58683641905a58f31c05605c953dc1d888466e8d434a6ceff47b76df03e"
        runtime_minutes: "1440"
    }
}

task CustomAnnotate {
    input {
        String pairName
        File vcfAnnotatedVep
        File cosmicCensus
        File cancerResistanceMutations
        String vcfAnnotatedCancerGeneCensusPath = "~{pairName}.v7.cosmic_census.vep.annotated.vcf"
        String vcfAnnotatedResistancePath = "~{pairName}.v7.resistance.vep.annotated.vcf"
        String vcfAnnotatedIdPath = "~{pairName}.v7.id.vep.annotated.vcf"
        String vcfCsqRenamedPath = "~{pairName}.snv.indel.supplemental.v7.CsqRenamedPath.vcf"
        String germFilteredOutFilePath = "~{pairName}.snv.indel.supplemental.v7.annotated.vcf"
        Int memoryGb = 2
        Int diskSize = ceil( size(cosmicCensus, "GB") * 2) + ceil( size(cancerResistanceMutations, "GB")) + ceil( size(vcfAnnotatedVep, "GB") * 2) + 5
        
    }

    command {
    set -e -o pipefail
    
    # AddCosmic
        python3 \
        /add_cancer_gene_census.py \
        ~{cosmicCensus} \
        ~{vcfAnnotatedVep} \
        ~{vcfAnnotatedCancerGeneCensusPath}
    # AddCancerResistanceMutations
        python3 \
        /add_cancer_resistance_mutations.py \
        ~{cancerResistanceMutations} \
        ~{vcfAnnotatedCancerGeneCensusPath} \
        ~{vcfAnnotatedResistancePath}
    # AnnotateId
        python3 \
        /annotate_id.py \
        ~{vcfAnnotatedResistancePath} \
        ~{vcfAnnotatedIdPath}
    # RenameCsqVcf
        python3 \
        /rename_csq_vcf.py \
        ~{vcfAnnotatedIdPath} \
        ~{vcfCsqRenamedPath}
    # GermFilterVcf
        python3 \
        /germ_filter_vcf.py \
        ~{vcfCsqRenamedPath} \
        ~{germFilteredOutFilePath}
    }
    
    output {
        File vcfCsqRenamed = "~{vcfCsqRenamedPath}"
    }
    
    runtime {
        mem: memoryGb + "G"
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-comp-s-fd4e/gcs_htslib_suite@sha256:47eba58683641905a58f31c05605c953dc1d888466e8d434a6ceff47b76df03e"
        runtime_minutes: "1440"
    }
}


task RemoveSpanning {
    input {
        String sampleId
        String noSpanningVcfPath = "~{sampleId}.v7.remove_spanning.vep.annotated.vcf"
        File vcfAnnotatedVep
        Int threads = 4
        Int memoryGb = 16
        Int diskSize = ceil( size(vcfAnnotatedVep, "GB") * 2) + 10
    }

    command {
        bcftools \
        view \
        --exclude 'ALT="*"' \
        --threads ~{threads} \
        -o ~{noSpanningVcfPath} \
        ~{vcfAnnotatedVep}
    }

    output {
        File noSpanningVcf = "~{noSpanningVcfPath}"
    }

    runtime {
        mem: memoryGb + "G"
        cpus: threads
        cpu : threads
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/bcftools@sha256:d9b84254b8cc29fcae76b728a5a9a9a0ef0662ee52893d8d82446142876fb400"
        runtime_minutes: "500"
    }
}

task AddCosmic {
    input {
        String pairName
        String vcfAnnotatedCancerGeneCensusPath = "~{pairName}.v7.cosmic_census.vep.annotated.vcf"
        File vcfAnnotatedVep
        File cosmicCensus
        Int memoryGb = 2
        Int diskSize = ceil( size(vcfAnnotatedVep, "GB") * 2) + 5
    }

    command {
        python \
        /add_cancer_gene_census.py \
        ~{cosmicCensus} \
        ~{vcfAnnotatedVep} \
        ~{vcfAnnotatedCancerGeneCensusPath}
    }

    output {
        File vcfAnnotatedCancerGeneCensus = "~{vcfAnnotatedCancerGeneCensusPath}"
    }

    runtime {
        mem: memoryGb + "G"
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/somatic_dna_tools@sha256:f281e73ddf515f3a5db6e766ef12fb2331dafa2b27c88e426764dcd8b4a7e19b"
        runtime_minutes: "1440"
    }
}

task AddCancerResistanceMutations {
    input {
        String pairName
        String vcfAnnotatedResistancePath = "~{pairName}.v7.resistance.vep.annotated.vcf"
        File vcfAnnotatedCancerGeneCensus
        File cancerResistanceMutations
        Int memoryGb = 2
        Int diskSize = ceil( size(vcfAnnotatedCancerGeneCensus, "GB") * 2) + ceil( size(cancerResistanceMutations, "GB")) + 5
    }

    command {
        python \
        /add_cancer_resistance_mutations.py \
        ~{cancerResistanceMutations} \
        ~{vcfAnnotatedCancerGeneCensus} \
        ~{vcfAnnotatedResistancePath}
    }

    output {
        File vcfAnnotatedResistance = "~{vcfAnnotatedResistancePath}"
    }

    runtime {
        mem: memoryGb + "G"
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/somatic_dna_tools@sha256:f281e73ddf515f3a5db6e766ef12fb2331dafa2b27c88e426764dcd8b4a7e19b"
        runtime_minutes: "1440"
    }
}

task AddCancerResistanceMutationsFinal {
    input {
        String pairName
        String vcfAnnotatedResistancePath = "~{pairName}.v7.resistance.vep.annotated.vcf"
        File vcfAnnotatedCancerGeneCensus
        File cancerResistanceMutations
        Int memoryGb = 2
        Int diskSize = ceil( size(vcfAnnotatedCancerGeneCensus, "GB") * 2) + ceil( size(cancerResistanceMutations, "GB")) + 5
    }

    command {
        python \
        /add_cancer_resistance_mutations.py \
        ~{cancerResistanceMutations} \
        ~{vcfAnnotatedCancerGeneCensus} \
        ~{vcfAnnotatedResistancePath}
    }

    output {
        File vcfAnnotatedResistance = "~{vcfAnnotatedResistancePath}"
    }

    runtime {
        mem: memoryGb + "G"
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/somatic_dna_tools@sha256:f281e73ddf515f3a5db6e766ef12fb2331dafa2b27c88e426764dcd8b4a7e19b"
        runtime_minutes: "1440"
    }
}

task AnnotateId {
    input {
        Int memoryGb = 2
        String pairName
        String vcfAnnotatedIdPath = "~{pairName}.v7.id.vep.annotated.vcf"
        File vcfAnnotatedResistance
        Int diskSize = ceil( size(vcfAnnotatedResistance, "GB") * 2) + 5
    }

    command {
        python \
        /annotate_id.py \
        ~{vcfAnnotatedResistance} \
        ~{vcfAnnotatedIdPath}
    }

    output {
        File vcfAnnotatedId = "~{vcfAnnotatedIdPath}"
    }

    runtime {
        mem: memoryGb + "G"
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/somatic_dna_tools@sha256:f281e73ddf515f3a5db6e766ef12fb2331dafa2b27c88e426764dcd8b4a7e19b"
        runtime_minutes: "500"
    }
}

task RenameCsqVcf {
    input {
        Int memoryGb = 2
        String pairName
        String vcfCsqRenamedPath = "~{pairName}.snv.indel.supplemental.v7.annotated.vcf"
        File vcfAnnotatedId
        Int diskSize = ceil( size(vcfAnnotatedId, "GB") * 2) + 5
    }

    command {
        python \
        /rename_csq_vcf.py \
        ~{vcfAnnotatedId} \
        ~{vcfCsqRenamedPath}
    }

    output {
        File vcfCsqRenamed = "~{vcfCsqRenamedPath}"
    }

    runtime {
        mem: memoryGb + "G"
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/somatic_dna_tools@sha256:f281e73ddf515f3a5db6e766ef12fb2331dafa2b27c88e426764dcd8b4a7e19b"
        runtime_minutes: "500"
    }
}

task MainVcf {
    input {
        Int memoryGb = 2
        String pairName
        String mainVcfPath = "~{pairName}.snv.indel.final.v7.annotated.vcf"
        File vcfAnnotated
        Int diskSize = ceil( size(vcfAnnotated, "GB") * 2) + 5
    }

    command {
        python \
        /make_main_vcf.py \
        ~{vcfAnnotated} \
        ~{mainVcfPath}
    }

    output {
        File mainVcf = "~{mainVcfPath}"
    }

    runtime {
        mem: memoryGb + "G"
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/somatic_dna_tools@sha256:f281e73ddf515f3a5db6e766ef12fb2331dafa2b27c88e426764dcd8b4a7e19b"
        runtime_minutes: "60"
    }
}

task TableVcf {
    input {
        Int memoryGb = 2
        String pairName
        String vcfAnnotatedTxtPath = "~{pairName}.snv.indel.final.v7.annotated.txt"
        String tumor
        String normal
        File mainVcf
        Int diskSize = ceil( size(mainVcf, "GB") * 2) + 5
    }

    command {
        /make_txt.py \
        --vcf ~{mainVcf} \
        --txt ~{vcfAnnotatedTxtPath} \
        --tumor ~{tumor} \
        --normal ~{normal}
    }

    output {
        File vcfAnnotatedTxt = "~{vcfAnnotatedTxtPath}"
    }

    runtime {
        mem: memoryGb + "G"
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/somatic_dna_tools@sha256:f281e73ddf515f3a5db6e766ef12fb2331dafa2b27c88e426764dcd8b4a7e19b"
        runtime_minutes: "60"
    }
}

task VcfToMultiCenterMaf {
    input {
        Int memoryGb = 2
        String pairName
        String intersectionMafPath
        String unionMafPath
        String library
        String vepGenomeBuild
        String tumor
        String normal
        File mainVcf
        Int diskSize = ceil( size(mainVcf, "GB") * 2) + 5
#         File makeMaf = "pipelines/snv-indel-concordance/somatic_dna_tools/make_maf_multicenter.py"
        File makeMaf = "gs://nygc-comp-s-fd4e-input/somatic/hcmi/tools/make_maf_multicenter.py"
    }

    command {
        python \
        ~{makeMaf} \
        --vcf ~{mainVcf} \
        --maf ~{unionMafPath} \
        --intersection-maf ~{intersectionMafPath} \
        --library ~{library} \
        --vep-version ~{vepGenomeBuild} \
        --tumor ~{tumor} \
        --normal ~{normal}
    }

    output {
        File unionMaf = "~{unionMafPath}"
        File intersectionMaf = "~{intersectionMafPath}"
    }

    runtime {
        mem: memoryGb + "G"
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/somatic_dna_tools@sha256:f281e73ddf515f3a5db6e766ef12fb2331dafa2b27c88e426764dcd8b4a7e19b"
        runtime_minutes: "90"
    }
}


task VcfToMaf {
    input {
        Int memoryGb = 2
        String pairName
        String mafPath = "~{pairName}.snv.indel.final.v7.annotated.maf"
        String library
        String vepGenomeBuild
        String tumor
        String normal
        File ensemblEntrez
        File mainVcf
        Int diskSize = ceil( size(mainVcf, "GB") * 2) + 5
    }

    command {
        python \
        /make_maf.py \
        --vcf ~{mainVcf} \
        --maf ~{mafPath} \
        --library ~{library} \
        --vep-version ~{vepGenomeBuild} \
        --tumor ~{tumor} \
        --normal ~{normal} \
        --ensembl-entrez ~{ensemblEntrez}
    }

    output {
        File maf = "~{mafPath}"
    }

    runtime {
        mem: memoryGb + "G"
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/somatic_dna_tools@sha256:f281e73ddf515f3a5db6e766ef12fb2331dafa2b27c88e426764dcd8b4a7e19b"
        runtime_minutes: "90"
    }
}

task annotateBicSeq2Cnv {
    input {
        Int memoryGb = 36
        String pairName
        Array[String] listOfChroms
        String tumor
        String normal
        File bicseq2
        File cytoBand
        File dgv
        File thousandG
        File cosmicUniqueBed
        File cancerCensusBed
        File ensemblUniqueBed

        String cnvAnnotatedFinalBedPath  = "~{pairName}.cnv.annotated.v7.final.bed"
        String cnvAnnotatedSupplementalBedPath  = "~{pairName}.cnv.annotated.v7.supplemental.bed"
        Int diskSize = 4
    }

    command {
        Rscript \
        /annotate-cnv.r \
        --cnv=~{bicseq2} \
        --caller=bicseq2 \
        --tumor=~{tumor} \
        --normal=~{normal} \
        --cytoband=~{cytoBand} \
        --db_names=DGV,1000G,COSMIC \
        --db_files=~{dgv},~{thousandG},~{cosmicUniqueBed} \
        --cancer_census=~{cancerCensusBed} \
        --ensembl=~{ensemblUniqueBed} \
        --allowed_chr=~{sep="," listOfChroms} \
        --overlap_fraction=0.8 \
        --out_file_main=~{cnvAnnotatedFinalBedPath} \
        --out_file_supplemental=~{cnvAnnotatedSupplementalBedPath}
    }

    output {
        File cnvAnnotatedFinalBed = "~{cnvAnnotatedFinalBedPath}"
        File cnvAnnotatedSupplementalBed = "~{cnvAnnotatedSupplementalBedPath}"
    }

    runtime {
        mem: memoryGb + "G"
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/sv_cnv@sha256:1c14a50d131323a2a4bab323cf224879776af8de37f93df79292fd2e63269274"
        runtime_minutes: "60"
    }
}

task mergeSv {
    input {
        Int memoryGb = 8
        String pairName
        Array[String] listOfChroms
        String tumor
        String normal
        File filteredMantaSV
        IndexedVcf gridssVcf
        String vepGenomeBuild
        Int slop = 300
        String svMergedFinalBedPePath  = "~{pairName}.sv.merged.v7.final.bedpe"
        String svMergedSupplementalBedPePath  = "~{pairName}.sv.merged.v7.supplemental.bedpe"
        Int diskSize = 4
        Int minSvLength = 500

    }

    command {
        Rscript \
        /merge-caller-vcfs.r \
        --vcf=~{filteredMantaSV},~{gridssVcf.vcf} \
        --caller=manta,gridss \
        --tumor=~{tumor} \
        --normal=~{normal} \
        --build=~{vepGenomeBuild} \
        --slop=~{slop} \
        --allowed_chr=~{sep="," listOfChroms} \
        --min_sv_length=~{minSvLength} \
        --out_file=~{svMergedFinalBedPePath} \
        --out_file_supplemental=~{svMergedSupplementalBedPePath}
    }

    output {
        File svMergedFinalBedPe = "~{svMergedFinalBedPePath}"
        File svMergedSupplementalBedPe = "~{svMergedSupplementalBedPePath}"
    }

    runtime {
        mem: memoryGb + "G"
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/sv_cnv@sha256:1c14a50d131323a2a4bab323cf224879776af8de37f93df79292fd2e63269274"
        runtime_minutes: "180"
    }
}

task annotateSv {
    input {
        Int memoryGb = 8
        String pairName
        String tumor
        String normal

        Int slop = 500

        # gap,DGV,1000G,PON,COSMIC
        File gap
        File dgvBedpe
        File thousandGVcf
        File svPon
        File cosmicBedPe

        File svMergedBedPe
        String svMergedAnnotatedBedPePath  = "~{pairName}.sv.merged.annotated.v7.bedpe"
        Int diskSize = 4

    }

    command {
        Rscript \
        /annotate-bedpe-with-databases.r \
        --db_names=gap,DGV,1000G,PON,COSMIC \
        --db_files=~{gap},~{dgvBedpe},~{thousandGVcf},~{svPon},~{cosmicBedPe} \
        --slop=~{slop} \
        --db_ignore_strand=COSMIC \
        --bedpe=~{svMergedBedPe} \
        --out_file=~{svMergedAnnotatedBedPePath}
    }

    output {
        File svMergedAnnotatedBedPe = "~{svMergedAnnotatedBedPePath}"
    }

    runtime {
        mem: memoryGb + "G"
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/sv_cnv@sha256:1c14a50d131323a2a4bab323cf224879776af8de37f93df79292fd2e63269274"
        runtime_minutes: "180"
    }
}

task annotateGenesSv {
    input {
        Int memoryGb = 8
        String pairName
        String tumor
        String normal

        File ensemblUniqueBed
        File cancerCensusBed

        File svMergedAnnotatedFinalBedPe
        String svGeneAnnotatedFinalBedPePath  = "~{pairName}.sv.merged.gene.annotated.v7.bedpe"
        Int diskSize = 4

    }

    command {
        Rscript \
        /annotate-bedpe-with-genes.r \
        --ensembl=~{ensemblUniqueBed} \
        --cancer_census=~{cancerCensusBed} \
        --bedpe=~{svMergedAnnotatedFinalBedPe} \
        --out_file=~{svGeneAnnotatedFinalBedPePath}
    }

    output {
        File svGeneAnnotatedFinalBedPe = "~{svGeneAnnotatedFinalBedPePath}"
    }

    runtime {
        mem: memoryGb + "G"
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/sv_cnv@sha256:1c14a50d131323a2a4bab323cf224879776af8de37f93df79292fd2e63269274"
        runtime_minutes: "1440"
    }
}

task annotateGenesSvSupplemental {
    input {
        Int memoryGb = 8
        String pairName
        String tumor
        String normal

        File ensemblUniqueBed
        File cancerCensusBed

        File svMergedAnnotatedSupplementalBedPe
        String svGeneAnnotatedSupplementalBedPePath  = "~{pairName}.sv.merged.gene.annotated.v7.supplemental.bedpe"
        Int diskSize = 4
    }

    command {
        Rscript \
        /annotate-bedpe-with-genes.r \
        --ensembl=~{ensemblUniqueBed} \
        --cancer_census=~{cancerCensusBed} \
        --bedpe=~{svMergedAnnotatedSupplementalBedPe} \
        --out_file=~{svGeneAnnotatedSupplementalBedPePath} \
        --supplemental
    }

    output {
        File svGeneAnnotatedSupplementalBedPe = "~{svGeneAnnotatedSupplementalBedPePath}"
    }

    runtime {
        mem: memoryGb + "G"
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/sv_cnv@sha256:1c14a50d131323a2a4bab323cf224879776af8de37f93df79292fd2e63269274"
        runtime_minutes: "1440"
    }
}

task annotateWithCnvSv {
    input {
        Int memoryGb = 8
        String pairName
        String tumor
        String normal

        File cnvAnnotatedFinalBed

        File svGeneAnnotatedBedPe
        String svCnvAnnotatedBedPePath = "~{pairName}.sv.merged.gene.cnv.annotated.v7.bedpe"
        Int diskSize = 4
    }

    command {
        Rscript \
        /annotate-bedpe-with-cnv.r \
        --cnv=~{cnvAnnotatedFinalBed} \
        --bedpe=~{svGeneAnnotatedBedPe} \
        --out_file=~{svCnvAnnotatedBedPePath}
    }

    output {
        File svCnvAnnotatedBedPe = "~{svCnvAnnotatedBedPePath}"
    }

    runtime {
        mem: memoryGb + "G"
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/sv_cnv@sha256:1c14a50d131323a2a4bab323cf224879776af8de37f93df79292fd2e63269274"
        runtime_minutes: "180"
    }
}

task filterBedPe {
    input {
        Int memoryGb = 8
        String pairName

        File svCnvAnnotatedBedPe
        String svBedPePath
        String svHighConfidenceBedPePath
        Int diskSize = 4

    }

    command {
        Rscript \
        /filter-bedpe.r \
        --max_changepoint_distance=1000 \
        --filter_databases=DGV,1000G,PON \
        --bedpe=~{svCnvAnnotatedBedPe} \
        --out_file_somatic=~{svBedPePath} \
        --out_file_highconf=~{svHighConfidenceBedPePath}
    }

    output {
        File svBedPe = "~{svBedPePath}"
        File svHighConfidenceBedPe = "~{svHighConfidenceBedPePath}"
    }

    runtime {
        mem: memoryGb + "G"
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/sv_cnv@sha256:1c14a50d131323a2a4bab323cf224879776af8de37f93df79292fd2e63269274"
        runtime_minutes: "180"
    }
}
