version 1.0

import "../wdl_structs.wdl"

# General tasks

task PrepNygcVcf    {
    input {
        File nygcVcf
        Array[String] normals
        String nygcHighConfVcfPath = sub(basename(nygcVcf), ".vcf$", ".highConfidence.vcf")
        String nygcPreppedVcfPath = sub(basename(nygcVcf), ".vcf$", ".prepped.vcf")
        #File addGT = "gs://nygc-comp-s-fd4e-input/hcmi/tools/add_gt.py"
        File addGT = "pipelines/snv-indel-concordance/somatic_dna_tools/add_gt.py"
        Int? diskSize = (ceil( size(nygcVcf, "GB") )  * 2 ) + 4
        Int? memoryGb = 8
    }
    Int jvmHeap = memoryGb * 750  # Heap size in Megabytes. mem is in GB. (75% of mem)

    command {
        cat ~{nygcVcf} \
        | grep "^#" \
        > ~{nygcHighConfVcfPath}

        cat ~{nygcVcf} \
        | grep -v "^#" \
        | grep "HighConfidence" \
        >> ~{nygcHighConfVcfPath}

        python \
        ~{addGT} \
        --vcf ~{nygcHighConfVcfPath} \
        --normals ~{sep=" " normals} \
        --out ~{nygcPreppedVcfPath}

    }

    output {
        File nygcPreppedVcf = "~{nygcPreppedVcfPath}"
    }

    runtime {
        mem: memoryGb + "G"
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/somatic_dna_tools@sha256:f281e73ddf515f3a5db6e766ef12fb2331dafa2b27c88e426764dcd8b4a7e19b"
        runtime_minutes: "720"
    }
}

task PrepWashuVcf    {
    input {
        File washuVcf
        String normal
        String tumor
        String washuPreppedVcfPath = sub(basename(washuVcf), ".vcf$", ".prepped.vcf")
        # File prepWashu = "gs://nygc-comp-s-fd4e-input/hcmi/tools/prep_washu.py"
        File prepWashu = "pipelines/snv-indel-concordance/somatic_dna_tools/prep_washu.py"
        Int? diskSize = (ceil( size(washuVcf, "GB") )  * 2 ) + 4
        Int? memoryGb = 8
    }
    command {
        python ${prepWashu} \
        ${washuVcf} \
        ${washuPreppedVcfPath} \
        ${normal} \
        ${tumor}
    }

    output {
        File washuPreppedVcf = "~{washuPreppedVcfPath}"
    }

    runtime {
        mem: memoryGb + "G"
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/somatic_dna_tools@sha256:f281e73ddf515f3a5db6e766ef12fb2331dafa2b27c88e426764dcd8b4a7e19b"
        runtime_minutes: "720"
    }
}

task PrepBroadMaf    {
    input {
        File broadMaf
        File referenceFa
        String normal
        String tumor
        String cleanMaf = sub(basename(broadMaf), "maf$", ".cleaned.maf")
        String cleanVcf = sub(basename(broadMaf), "maf$", ".cleaned.vcf")
        # File convertMaf = "gs://nygc-comp-s-fd4e-input/hcmi/tools/maf2vcf.pl"
        File convertMaf = "pipelines/snv-indel-concordance/somatic_dna_tools/maf2vcf.pl"
        Int? diskSize = (ceil( size(broadMaf, "GB") )  * 2 ) + 4
        Int? memoryGb = 8
    }
    command {
        cat ${broadMaf} \
        | sed 's|__UNKNOWN__||g' \
        > ${cleanMaf}

        mkdir out/

        perl ${convertMaf} \
        --input-maf ${cleanMaf} \
        --output-dir out/ \
        --ref-fasta ${referenceFa}
    }

    output {
        File broadCleanVcf = "out/~{cleanVcf}"
    }

    runtime {
        mem: memoryGb + "G"
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/samtools@sha256:963b0b2f24908832efab8ddccb7a7f3ba5dca9803bc099be7cf3a455766610fd"
        runtime_minutes: "720"
    }
}

task RenameBroadVcf    {
    input {
        File cleanVcf
        String normal
        String tumor
        String renamedVcf = sub(basename(cleanVcf), ".cleaned.vcf$", ".renamed.vcf")
        #File prepBroad = "gs://nygc-comp-s-fd4e-input/hcmi/tools/prep_broad.py"
        File prepBroad = "pipelines/snv-indel-concordance/somatic_dna_tools/prep_broad.py"
        Int? diskSize = (ceil( size(cleanVcf, "GB") )  * 2 ) + 4
        Int? memoryGb = 8
    }
    command {
        python ${prepBroad} \
        ${cleanVcf} \
        ${renamedVcf} \
        ${normal} \
        ${tumor}
    }

    output {
        File broadVcf = "~{renamedVcf}"
    }

    runtime {
        mem: memoryGb + "G"
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/somatic_dna_tools@sha256:f281e73ddf515f3a5db6e766ef12fb2331dafa2b27c88e426764dcd8b4a7e19b"
        runtime_minutes: "720"
    }
}

task PrepExomeWgsMetadata {
    input {
        String pairName
        File callerVcf
        String renameMetaVcfPath = "~{pairName}.rename_metadata.vcf"
        String renameVcfPath = "~{pairName}.rename.vcf"
        String prepCallerVcfPath = "~{pairName}.merge_prep.vcf"
        String splitVcfPath = "~{pairName}.split.vcf"
        String mnvVcfPath = "~{pairName}.mnv.vcf"
        IndexedReference referenceFa
        String tool
        String normalId
        String tumorId
        Int memoryGb = 16
        Int threads = 4
        Int diskSize = (ceil( size(callerVcf, "GB") )  * 2 ) + 4
    }

    command {
        set -e -o pipefail
        # RenameExomeWgsMetadata
        python \
        /rename_exome_wgs_metadata.py \
        ~{callerVcf} \
        ~{renameMetaVcfPath} \
        ~{tool}
        # MergePrepExomeWgs
        python \
        /merge_exome_wgs_prep.py \
        --vcf ~{renameMetaVcfPath} \
        --out ~{prepCallerVcfPath} \
        --tool ~{tool}
        # RenameVcf
        python \
        /rename_vcf.py \
        ~{prepCallerVcfPath} \
        ~{renameVcfPath} \
        ~{normalId} \
        ~{tumorId} \
        ~{tool}
        
#        CompressIndexVcf
        bgzip -c \
        ~{renameVcfPath} \
        > ~{renameVcfPath}.gz
        
        tabix \
        -p vcf \
        ~{renameVcfPath}.gz
        
        # task SplitMultiAllelic
        bcftools \
        norm \
        -m \
        -any \
        --threads ~{threads} \
        --no-version \
        -f ~{referenceFa.fasta} \
        -o ~{splitVcfPath} \
        ~{renameVcfPath}.gz
    # task SplitMnv {
        python3 \
        /split_mnv.py \
        ~{splitVcfPath} \
        ~{mnvVcfPath} \
        ~{tool}
    }

    output {
        File mnvVcf = "~{mnvVcfPath}"
    }

    runtime {
        mem: memoryGb + "G"
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-comp-s-fd4e/gcs_htslib_suite@sha256:47eba58683641905a58f31c05605c953dc1d888466e8d434a6ceff47b76df03e"
        runtime_minutes: "720"
    }
}

task FilterPairHighConfidence  {
    input {
        Int memoryGb = 4
        String outputVcfPath
        File intersectionPairVcf
        IndexedVcf intersectionParticipantVcf
        String tumorBarcodeAliquot
        String normalBarcodeAliquot
        String participantId
        
        Int diskSize
        File makePairVcfs = "pipelines/snv-indel-concordance/somatic_dna_tools/make_pair_vcfs.py"
#         File makePairVcfs = "gs://nygc-comp-s-fd4e-input/hcmi/tools/make_pair_vcfs.py"
    }

    command {
        set -e -o pipefail


        python ~{makePairVcfs} \
        --normal ~{normalBarcodeAliquot} \
        --participant-vcf ~{intersectionParticipantVcf.vcf} \
        --pair-vcf ~{intersectionPairVcf} \
        --out-vcf ~{outputVcfPath}
    }

    output {
        File intersectionPairHighConfidencVcf = "~{outputVcfPath}"
    }

    runtime {
        mem: memoryGb + "G"
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-comp-s-fd4e/gcs_htslib_suite@sha256:47eba58683641905a58f31c05605c953dc1d888466e8d434a6ceff47b76df03e"
        runtime_minutes: "1440"
    }
}

task MakeUnionVcf {
    input {
        Int memoryGb = 2
        String pairName
        String intersectionVcf
        String unionVcfPath = "~{pairName}.snv.indel.union.vcf"
        Int diskSize = ceil( size(intersectionVcf, "GB") * 2) + 5
        File multiPatientAnnotation = "pipelines/snv-indel-concordance/somatic_dna_tools/multi_patient_annotation.py"
#         File multiPatientAnnotation = "gs://nygc-comp-s-fd4e-input/somatic/hcmi/tools/multi_patient_annotation.py"
    }

    command {
        
        python \
        ~{multiPatientAnnotation} \
        --vcf ~{intersectionVcf} \
        --multi-caller-out ~{unionVcfPath}
    }

    output {
        File unionVcf = "~{unionVcfPath}"
    }

    runtime {
        mem: memoryGb + "G"
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-comp-s-fd4e/gcs_htslib_suite@sha256:47eba58683641905a58f31c05605c953dc1d888466e8d434a6ceff47b76df03e"
        runtime_minutes: "1440"
    }
}

task MergeBcftools {
    input {
        Array[File] vcfs
        String mergedVcfPath = "mergedVcfPath.vcf"
        String sortedVcfPath
        String dollarSign = "$"
        Int threads = 1
        Int memoryGb = 16
        Int diskSize = (ceil( size(vcfs, "GB") )  * 3 ) + 10
    }

    command {
        grep \
         "^#" \
        ~{vcfs[0]} \
        > ~{mergedVcfPath}
        
        for vcf_file in ~{ sep=' ' vcfs } ; do
            grep \
            -v "^#" \
            ~{dollarSign}vcf_file \
            >> ~{mergedVcfPath}
        done
        
        bgzip -c \
        ~{mergedVcfPath} \
        > ~{mergedVcfPath}.gz
        
        tabix \
        -p vcf \
        ~{mergedVcfPath}.gz
        
        bcftools \
        sort \
        -o ~{sortedVcfPath} \
        ~{mergedVcfPath}.gz
    }

    output {
        File sortedVcf = "~{sortedVcfPath}"
    }

    runtime {
        mem: memoryGb + "G"
        cpus: threads
        cpu : threads
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-comp-s-fd4e/gcs_htslib_suite@sha256:47eba58683641905a58f31c05605c953dc1d888466e8d434a6ceff47b76df03e"
        runtime_minutes: "90"
    }
}

task CountAliquots {
    input {
        File mergedVcf
        String participantId
        String splitId
        Bam normalFinalBam
        Array[Bam] tumorFinalBams
        Array[File] tumorFinalBamList
        # intermediate files
        String countedVcfPath = "~{participantId}.splitId.counted.vcf"
        String countedOrderedVcfPath = "~{participantId}.splitId.counted.ordered.vcf"
        Int threads = 1
        Int memoryGb = 16
        Int diskSize = (ceil( size(normalFinalBam.bam, "GB") + size(tumorFinalBamList, "GB"))) + 4
        # public copy avilable in
        # https://bitbucket.nygenome.org/projects/WDL/repos/somatic_dna_wdl/browse/tools/rename_center_metadata.py?at=refs%2Fheads%2Fcenter_merge
        File addNygcAlleleCounts = "pipelines/snv-indel-concordance/somatic_dna_tools/add_nygc_allele_counts_to_vcf_participant.py"
#         File addNygcAlleleCounts = "gs://nygc-comp-s-fd4e-input/somatic/hcmi/tools/add_nygc_allele_counts_to_vcf_participant.py"
        File reOrderGt = "pipelines/snv-indel-concordance/somatic_dna_tools/re_order_gt.py"
#         File reOrderGt = "gs://nygc-comp-s-fd4e-input/somatic/hcmi/tools/re_order_gt.py"

    }

    command {
    set -e -o pipefail
    
    # task addNygcAlleleCounts
    python ~{addNygcAlleleCounts} \
    --vcf ~{mergedVcf} \
    --tumor-bams ~{ sep=' ' tumorFinalBamList } \
    --normal-bam ~{normalFinalBam.bam} \
    --out ~{countedVcfPath}
    
    # task reOrderGt
    # avoid occasional bug where GATK breaks lines if no GT and 
    # format item that ends in GT is not first 
    python \
    ~{reOrderGt} \
    ~{countedVcfPath} \
    ~{countedOrderedVcfPath}
    
    }

    output {
        File countedVcf = "~{countedOrderedVcfPath}"
    }

    runtime {
        mem: memoryGb + "G"
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-comp-s-fd4e/gcs_htslib_suite@sha256:47eba58683641905a58f31c05605c953dc1d888466e8d434a6ceff47b76df03e"
        runtime_minutes: "720"
    }
}


task AnnotateWithMutect2 {
    input {
        File mergedVcf
        File mutect2
        String participantId
        File cosmicCensus
        # intermediate files
        String mutect2AnnotatedVcfPath = "~{participantId}.merged.mutect2_annotated.vcf"
        String mutect2AnnoConsolidatedVcfPath = "~{participantId}.mutect2_annotated_consolidated.vcf"
        String mutect2AnnoConsolidatedMnvPath = "~{participantId}.mutect2_annotated_mnv.vcf"
        String vcfAnnotatedCancerGeneCensusPath = "~{participantId}.cosmic.vcf"
        String intersectionVcfPath = "~{participantId}.intersection.vcf"
        String unionPath = "~{participantId}.union.vcf"
        Int threads = 4
        Int memoryGb = 16
        Int diskSize = (ceil( size(mergedVcf, "GB") )  * 4 ) + 50
        # public copy avilable in
        # https://bitbucket.nygenome.org/projects/WDL/repos/somatic_dna_wdl/browse/tools/rename_center_metadata.py?at=refs%2Fheads%2Fcenter_merge
        File mutect2Annotation = "pipelines/snv-indel-concordance/somatic_dna_tools/mutect2_annotation.py"
        File consolidateCsq = "pipelines/snv-indel-concordance/somatic_dna_tools/consolidate_csq.py"
        File consolidateMnvs = "pipelines/snv-indel-concordance/somatic_dna_tools/consolidate_mnvs.py"
        File multiPatientAnnotation = "pipelines/snv-indel-concordance/somatic_dna_tools/multi_patient_annotation.py"
        File reOrderGt = "pipelines/snv-indel-concordance/somatic_dna_tools/re_order_gt.py"
        File addCancerGeneCensus = "pipelines/snv-indel-concordance/somatic_dna_tools/add_cancer_gene_census.py"

#         File reOrderGt = "gs://nygc-comp-s-fd4e-input/somatic/hcmi/tools/re_order_gt.py"
#         File multiPatientAnnotation = "gs://nygc-comp-s-fd4e-input/somatic/hcmi/tools/multi_patient_annotation.py"
#         File mutect2Annotation = "gs://nygc-comp-s-fd4e-input/somatic/hcmi/tools/mutect2_annotation.py"
#         File consolidateCsq = "gs://nygc-comp-s-fd4e-input/somatic/hcmi/tools/consolidate_csq.py"
#         File consolidateMnvs = "gs://nygc-comp-s-fd4e-input/somatic/hcmi/tools/consolidate_mnvs.py"
#         File addCancerGeneCensus = "gs://nygc-comp-s-fd4e-input/somatic/hcmi/tools/add_cancer_gene_census.py"

    }

    command {
    set -e -o pipefail
    
    # task IndexCompressVcf
        bgzip -c \
        ~{mutect2} \
        > ~{mutect2}.gz
        tabix \
        -p vcf \
        ~{mutect2}.gz
    # task Mutect2Annotation
    python ~{mutect2Annotation} \
    --vcf ~{mergedVcf} \
    --mutect2-vcf ~{mutect2}.gz \
    --out ~{mutect2AnnotatedVcfPath}
    
    # task consolidateCsq
    python ~{consolidateCsq} \
    --vcf ~{mutect2AnnotatedVcfPath} \
    --out ~{mutect2AnnoConsolidatedVcfPath}
    
    
    # task consolidateMnvs
    python ~{consolidateMnvs} \
    --vcf ~{mutect2AnnoConsolidatedVcfPath} \
    --out ~{mutect2AnnoConsolidatedMnvPath}
    # AddCosmic
    python3 \
    ~{addCancerGeneCensus} \
    ~{cosmicCensus} \
    ~{mutect2AnnoConsolidatedMnvPath} \
    ~{vcfAnnotatedCancerGeneCensusPath}
    # task reOrderGt
    # avoid occasional bug where GATK breaks lines if no GT and 
    # format item that ends in GT is not first 
    python \
    ~{reOrderGt} \
    ~{vcfAnnotatedCancerGeneCensusPath} \
    ~{unionPath}
    # task MakeIntersectionVcf
    python \
    ~{multiPatientAnnotation} \
    --vcf ~{unionPath} \
    --multi-caller-out ~{intersectionVcfPath}
    }

    output {
        File unionVcf = "~{unionPath}"
        File intersectionVcf = "~{intersectionVcfPath}"
    }

    runtime {
        mem: memoryGb + "G"
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-comp-s-fd4e/gcs_htslib_suite@sha256:47eba58683641905a58f31c05605c953dc1d888466e8d434a6ceff47b76df03e"
        runtime_minutes: "720"
    }
}


task SelectVcfSamples {
    input {
        Array[String] retainedSamples
        String nonNormalVcfPath  = "retainedSamples.vcf"
        File centerVcf
        IndexedReference referenceFa
        Int threads = 1
        Int memoryGb = 16
        Int diskSize = (ceil( size(centerVcf, "GB") )  * 3 ) + 10
    }

    command {
        bcftools \
        view \
        --threads ~{threads} \
        --no-version \
        --samples ~{sep="," retainedSamples} \
        -f ~{referenceFa.fasta} \
        -o ~{nonNormalVcfPath} \
        ~{centerVcf}
    }

    output {
        File nonNormalVcf = "~{nonNormalVcfPath}"
    }

    runtime {
        mem: memoryGb + "G"
        cpus: threads
        cpu : threads
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/bcftools@sha256:d9b84254b8cc29fcae76b728a5a9a9a0ef0662ee52893d8d82446142876fb400"
        runtime_minutes: "90"
    }
}

task WipeAnnnotation {
    input {
        Int memoryGb = 2
        File inputVcf
        String outputVcfPath
        String tmpVcfHeaderPath = "vcf_header_lines.hdr"
        Int diskSize = ceil( size(inputVcf, "GB") * 2) + 5
    }

    command {
        set -e -o pipefail


        bcftools annotate \
        --remove "INFO/CSQ,INFO" \
        ~{inputVcf} \
        > ~{outputVcfPath}
    }

    output {
        File outputVcf = "~{outputVcfPath}"
    }

    runtime {
        mem: memoryGb + "G"
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-comp-s-fd4e/gcs_htslib_suite@sha256:47eba58683641905a58f31c05605c953dc1d888466e8d434a6ceff47b76df03e"
        runtime_minutes: "1440"
    }
}

task PrepParticipantCalls {
    input {
        String tumorBarcodeAliquot
        File centerVcf
        IndexedReference referenceFa
        # intermediate files
        String renameMetaVcfPath = sub(basename(centerVcf), "$", ".rename_metadata.vcf")
        String prepCallerVcfPath = sub(basename(centerVcf), "$", ".merge_prep.vcf")
        String splitVcfPath = sub(basename(centerVcf), "$", ".split.vcf")
        String mnvVcfPath = sub(basename(centerVcf), "$", ".mnv.vcf")
        String mnvSortedVcfPath = sub(basename(centerVcf), "$", ".mnv.sorted.vcf")
        String nonNormalKind
        Int index
        Int threads = 4
        Int memoryGb = 16
        Int diskSize = (ceil( size(centerVcf, "GB") )  * 8 ) + 10
        # public copy avilable in
        # https://bitbucket.nygenome.org/projects/WDL/repos/somatic_dna_wdl/browse/tools/rename_center_metadata.py?at=refs%2Fheads%2Fcenter_merge
#         File renameAliquotMetadata = "pipelines/snv-indel-concordance/somatic_dna_tools/rename_aliquot_metadata.py"
#         File mergeAliquotPrep = "pipelines/snv-indel-concordance/somatic_dna_tools/merge_aliquot_prep.py"
        File renameAliquotMetadata = "gs://nygc-comp-s-fd4e-input/somatic/hcmi/tools/rename_aliquot_metadata.py"
        File mergeAliquotPrep = "gs://nygc-comp-s-fd4e-input/somatic/hcmi/tools/merge_aliquot_prep.py"
    }

    command {
    set -e -o pipefail
    
    # task RenameAliquotMetadata
        python3 \
        ~{renameAliquotMetadata} \
        ~{centerVcf} \
        ~{renameMetaVcfPath} \
        Aliquot~{index} \
        ~{tumorBarcodeAliquot}
        
    # task MergePrepSupport
        python3 \
        ~{mergeAliquotPrep} \
        --vcf ~{renameMetaVcfPath} \
        --out ~{prepCallerVcfPath} \
        --aliquot ~{tumorBarcodeAliquot} \
        --non-normal-kind ~{nonNormalKind}
    # task confirmedIndexCompressVcf
        bgzip -c \
        ~{prepCallerVcfPath} \
        > ~{prepCallerVcfPath}.gz
        tabix \
        -p vcf \
        ~{prepCallerVcfPath}.gz
    # task SplitMultiAllelic
        bcftools \
        norm \
        -m \
        -any \
        --threads ~{threads} \
        --no-version \
        -f ~{referenceFa.fasta} \
        -o ~{splitVcfPath} \
        ~{prepCallerVcfPath}.gz
    # task SplitMnv {
        python3 \
        /split_mnv.py \
        ~{splitVcfPath} \
        ~{mnvVcfPath} \
        ~{tumorBarcodeAliquot}
    
    # task mnvIndexCompressVcf
        grep \
        "^#" \
        ~{mnvVcfPath} \
        > ~{mnvSortedVcfPath} 
        grep \
        -v "^#" \
        ~{mnvVcfPath} \
        | sort -k1,1 -k2,2n \
        >> ~{mnvSortedVcfPath} 
        bgzip -c \
        ~{mnvSortedVcfPath} \
        > ~{mnvSortedVcfPath}.gz
        tabix \
        -p vcf \
        ~{mnvSortedVcfPath}.gz
    }

    output {
        IndexedVcf mnvVcf = object {
            vcf : "~{mnvSortedVcfPath}.gz",
            index : "~{mnvSortedVcfPath}.gz.tbi"
        }
    }

    runtime {
        mem: memoryGb + "G"
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-comp-s-fd4e/gcs_htslib_suite@sha256:47eba58683641905a58f31c05605c953dc1d888466e8d434a6ceff47b76df03e"
        runtime_minutes: "720"
    }
}


task PrepCenterCallsExomeWgs {
    input {
        String pairName
        File callerVcf
        IndexedReference referenceFa
        # intermediate files
        String renameMetaVcfPath = sub(basename(callerVcf), "$", ".rename_metadata.vcf")
        String prepCallerVcfPath = sub(basename(callerVcf), "$", ".merge_prep.vcf")
        String renameVcfPath = sub(basename(callerVcf), "$", ".rename.vcf")
        String splitVcfPath = sub(basename(callerVcf), "$", ".split.vcf")
        String mnvVcfPath = sub(basename(callerVcf), "$", ".mnv.vcf")
        String tool
        String normalId
        String tumorId
        Int threads = 4
        Int memoryGb = 16
        Int diskSize = (ceil( size(callerVcf, "GB") )  * 2 ) + 4
        # public copy avilable in
        # https://bitbucket.nygenome.org/projects/WDL/repos/somatic_dna_wdl/browse/tools/rename_center_metadata.py?at=refs%2Fheads%2Fcenter_merge
        File renameCenterMetadata = "gs://nygc-comp-s-fd4e-input/somatic/hcmi/tools/rename_center_metadata.py"
        File mergeCenterPrep = "gs://nygc-comp-s-fd4e-input/somatic/hcmi/tools/merge_center_prep.py"
#         File renameCenterMetadata = "pipelines/snv-indel-concordance/somatic_dna_tools/rename_center_metadata.py"
#         File mergeCenterPrep = "pipelines/snv-indel-concordance/somatic_dna_tools/merge_center_prep.py"
    }

    command {
    set -e -o pipefail
    
    # task RenameCenterMetadata
        python3 \
        ~{renameCenterMetadata} \
        ~{callerVcf} \
        ~{renameMetaVcfPath} \
        ~{tool}
    # task MergePrepSupport
        python3 \
        ~{mergeCenterPrep} \
        --exome-wgs \
        --vcf ~{renameMetaVcfPath} \
        --out ~{prepCallerVcfPath} \
        --tool ~{tool}
    # task RenameVcf
        python3 \
        /rename_vcf.py \
        ~{prepCallerVcfPath} \
        ~{renameVcfPath} \
        ~{normalId} \
        ~{tumorId} \
        ~{tool}
    # task confirmedIndexCompressVcf
        bgzip -c \
        ~{renameVcfPath} \
        > ~{renameVcfPath}.gz
        tabix \
        -p vcf \
        ~{renameVcfPath}.gz
    # task SplitMultiAllelic
        bcftools \
        norm \
        -m \
        -any \
        --threads ~{threads} \
        --no-version \
        -f ~{referenceFa.fasta} \
        -o ~{splitVcfPath} \
        ~{renameVcfPath}.gz
    # task SplitMnv {
        python3 \
        /split_mnv.py \
        ~{splitVcfPath} \
        ~{mnvVcfPath} \
        ~{tool}
    }

    output {
        File mnvVcf = "~{mnvVcfPath}"
    }

    runtime {
        mem: memoryGb + "G"
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-comp-s-fd4e/gcs_htslib_suite@sha256:47eba58683641905a58f31c05605c953dc1d888466e8d434a6ceff47b76df03e"
        runtime_minutes: "720"
    }
}

task MergeExomeWgsCallersFull {
    input {
        String chrom
        String pairName
        String tumorId
        String normalId
        String mergedChromVcfPath = "~{pairName}.merged_supported.v7.~{chrom}.vcf"
        String columnChromVcfPath = "~{pairName}.single_column.v7.~{chrom}.vcf"
        String finalChromVcfPath = "~{pairName}.mnv.final.v7.filtered.~{chrom}.vcf"
        Array[IndexedVcf] allVcfCompressed
        Array[File] allVcfCompressedList
        Int threads = 4
        Int memoryGb = 16
        Int diskSize = 20
    }

    command {
        # MergeExomeWgsCallers
        bcftools \
        merge \
        -r ~{chrom} \
        --force-samples \
        --no-version \
        --threads ~{threads} \
        -f PASS \
        -F x \
        -m none \
        -o ~{mergedChromVcfPath} \
        -i called_by:join,num_callers:sum,MNV_ID:join \
        ~{sep=" " allVcfCompressedList}
        
        # MergeColumns
        python3 \
        /merge_columns.py \
        ~{mergedChromVcfPath} \
        ~{columnChromVcfPath} \
        ~{tumorId} \
        ~{normalId}
        # SnvstomnvsAnnotateExomeWgsCalled
        python \
        /SNVsToMNVs_AnnotateExomeWgsCalled.py \
        -i ~{columnChromVcfPath} \
        -o ~{finalChromVcfPath}
    }

    output {
        File columnChromVcf = "~{columnChromVcfPath}"
        File finalChromVcf = "~{finalChromVcfPath}"
    }

    runtime {
        mem: memoryGb + "G"
        cpus: threads
        cpu : threads
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-comp-s-fd4e/gcs_htslib_suite@sha256:47eba58683641905a58f31c05605c953dc1d888466e8d434a6ceff47b76df03e"
        runtime_minutes: "500"
    }
}

task PrepCenterCalls {
    input {
        String pairName
        File callerVcf
        IndexedReference referenceFa
        # intermediate files
        String renameMetaVcfPath = sub(basename(callerVcf), "$", ".rename_metadata.vcf")
        String prepCallerVcfPath = sub(basename(callerVcf), "$", ".merge_prep.vcf")
        String renameVcfPath = sub(basename(callerVcf), "$", ".rename.vcf")
        String splitVcfPath = sub(basename(callerVcf), "$", ".split.vcf")
        String mnvVcfPath = sub(basename(callerVcf), "$", ".mnv.vcf")
        String tool
        String normalId
        String tumorId
        Int threads = 4
        Int memoryGb = 16
        Int diskSize = (ceil( size(callerVcf, "GB") )  * 2 ) + 4
        # public copy avilable in
        # https://bitbucket.nygenome.org/projects/WDL/repos/somatic_dna_wdl/browse/tools/rename_center_metadata.py?at=refs%2Fheads%2Fcenter_merge
        File renameCenterMetadata = "gs://nygc-comp-s-fd4e-input/somatic/hcmi/tools/rename_center_metadata.py"
        File mergeCenterPrep = "gs://nygc-comp-s-fd4e-input/somatic/hcmi/tools/merge_center_prep.py"
#         File renameCenterMetadata = "pipelines/snv-indel-concordance/somatic_dna_tools/rename_center_metadata.py"
#         File mergeCenterPrep = "pipelines/snv-indel-concordance/somatic_dna_tools/merge_center_prep.py"
    }

    command {
    set -e -o pipefail
    
    # task RenameCenterMetadata
        python3 \
        ~{renameCenterMetadata} \
        ~{callerVcf} \
        ~{renameMetaVcfPath} \
        ~{tool}
    # task MergePrepSupport
        python3 \
        ~{mergeCenterPrep} \
        --vcf ~{renameMetaVcfPath} \
        --out ~{prepCallerVcfPath} \
        --tool ~{tool}
    # task RenameVcf
        python3 \
        /rename_vcf.py \
        ~{prepCallerVcfPath} \
        ~{renameVcfPath} \
        ~{normalId} \
        ~{tumorId} \
        ~{tool}
    # task confirmedIndexCompressVcf
        bgzip -c \
        ~{renameVcfPath} \
        > ~{renameVcfPath}.gz
        tabix \
        -p vcf \
        ~{renameVcfPath}.gz
    # task SplitMultiAllelic
        bcftools \
        norm \
        -m \
        -any \
        --threads ~{threads} \
        --no-version \
        -f ~{referenceFa.fasta} \
        -o ~{splitVcfPath} \
        ~{renameVcfPath}.gz
    # task SplitMnv {
        python3 \
        /split_mnv.py \
        ~{splitVcfPath} \
        ~{mnvVcfPath} \
        ~{tool}
    }

    output {
        File mnvVcf = "~{mnvVcfPath}"
    }

    runtime {
        mem: memoryGb + "G"
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-comp-s-fd4e/gcs_htslib_suite@sha256:47eba58683641905a58f31c05605c953dc1d888466e8d434a6ceff47b76df03e"
        runtime_minutes: "720"
    }
}


task CompressIndexVcf {
    input {
        File vcf
        String vcfCompressedPath = sub(basename(vcf), ".vcf$", ".vcf.gz")
        String vcfIndexedPath = sub(basename(vcf), ".vcf$", ".vcf.gz.tbi")
        Int memoryGb = 8
        Int diskSize = (ceil( size(vcf, "GB") )  * 2 ) + 4
    }
    
    command {
        set -e -o pipefail
        
        bgzip -c \
        ~{vcf} \
        > ~{vcfCompressedPath}
        
        tabix \
        -p vcf \
        ~{vcfCompressedPath}
    }
    
    output {
        IndexedVcf vcfCompressedIndexed = object {
            vcf : "~{vcfCompressedPath}",
            index : "~{vcfIndexedPath}"
        }
    }
    
    runtime {
        mem: memoryGb + "G"
        memory : memoryGb + "GB"
        disks: "local-disk " + diskSize + " HDD"
        docker: "gcr.io/nygc-comp-s-fd4e/gcs_htslib_suite@sha256:47eba58683641905a58f31c05605c953dc1d888466e8d434a6ceff47b76df03e"
        runtime_minutes: "1440"
    }
    
    meta {
        internalOnly : "False, produces output for external or internal BED files"
    }
}

task SplitBed {
    input {
        Int minSplits = 40
        File bed
        String splitBedPrefixPath
        String additionalSuffix = ".bed"
        Array[String] splitBedPath
        Int memoryGb = 16
        Int diskSize = 20
    }
    
    # populate output paths with (how to handle additional files?)
#    Array[String] suffixes = ["00", "01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29", "30", "31", "32", "33", "34", "35", "36", "37", "38", "39", "40", "41", "42", "43", "44", "45", "46", "47", "48", "49", "50"]
#    scatter (index in range(length(suffixes))) {
#        String suffix = suffixes[index]
#        String splitBedPath = "~{splitBedPrefixPath}~{suffix}~{additionalSuffix}"
#    }
    
    command {
        set -e -o pipefail
        
        split \
        --number=l/~{minSplits} \
        -d \
        -a2 \
        --additional-suffix=~{additionalSuffix} \
        ~{bed} \
        ~{splitBedPrefixPath}
    }
    
    output {
        Array[File?] outputBeds = splitBedPath
    }
    
    runtime {
        mem: memoryGb + "G"
        memory : memoryGb + "GB"
        disks: "local-disk " + diskSize + " HDD"
        docker: "gcr.io/nygc-comp-s-fd4e/gcs_htslib_suite@sha256:47eba58683641905a58f31c05605c953dc1d888466e8d434a6ceff47b76df03e"
        runtime_minutes: "90"
    }
    
    meta {
        internalOnly : "False, produces output for external or internal BED files"
    }
}

task getIntersection {
    input {
        String pairName
        String intersectedVcfPath = "~{pairName}.wgs_in_exonic.vcf"
        File intervalListBed
        File wgsVcf
        Int memoryGb = 16
        Int diskSize = 20
    }

    command {
        bedtools \
        intersect \
        -header \
        -a ~{wgsVcf} \
        -b ~{intervalListBed} \
        -F 1 \
        > ~{intersectedVcfPath}
    }

    output {
        File intersectedVcf = "~{intersectedVcfPath}"
    }

    runtime {
        mem: memoryGb + "G"
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/bedtools@sha256:9e737f5c96c00cf3b813d419d7a7b474c4013c9aa9dfe704eb36417570c6474e"
        runtime_minutes: "90"
    }
}

task RenameExomeWgsMetadata {
    input {
        String pairName
        File callerVcf
        String renameMetaVcfPath = sub(basename(callerVcf), "$", ".rename_metadata.vcf")
        String tool
        Int memoryGb = 16
        Int diskSize = (ceil( size(callerVcf, "GB") )  * 2 ) + 4
    }

    command {
        python \
        /rename_exome_wgs_metadata.py \
        ~{callerVcf} \
        ~{renameMetaVcfPath} \
        ~{tool}
    }

    output {
        File renameMetaVcf = "~{renameMetaVcfPath}"
    }

    runtime {
        mem: memoryGb + "G"
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/somatic_dna_tools@sha256:f281e73ddf515f3a5db6e766ef12fb2331dafa2b27c88e426764dcd8b4a7e19b"
        runtime_minutes: "90"
    }
}

task MergePrepExomeWgs {
    input {
        String pairName
        File renameMetaVcf
        String prepCallerVcfPath = sub(basename(renameMetaVcf, ".gz"), ".rename_metadata.vcf$", ".merge_prep.vcf")
        String tool
        Int memoryGb = 16
        Int diskSize = (ceil( size(renameMetaVcf, "GB") )  * 2 ) + 4
    }

    command {
        python \
        /merge_exome_wgs_prep.py \
        --vcf ~{renameMetaVcf} \
        --out ~{prepCallerVcfPath} \
        --tool ~{tool}
    }

    output {
        File prepCallerVcf = "~{prepCallerVcfPath}"
    }

    runtime {
        mem: memoryGb + "G"
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/somatic_dna_tools@sha256:f281e73ddf515f3a5db6e766ef12fb2331dafa2b27c88e426764dcd8b4a7e19b"
        runtime_minutes: "90"
    }
}

task MergeExomeWgsCallers {
    input {
        String chrom
        String pairName
        String mergedChromVcfPath = "~{pairName}.merged_supported.v7.~{chrom}.vcf"
        Array[IndexedVcf] allVcfCompressed
        Array[File] allVcfCompressedList
        Int threads = 8
        Int memoryGb = 16
        Int diskSize = 20
    }

    command {
        bcftools \
        merge \
        -r ~{chrom} \
        --force-samples \
        --no-version \
        --threads ~{threads} \
        -f PASS \
        -F x \
        -m none \
        -o ~{mergedChromVcfPath} \
        -i called_by:join,num_callers:sum,MNV_ID:join \
        ~{sep=" " allVcfCompressedList}
    }

    output {
        File mergedChromVcf = "~{mergedChromVcfPath}"
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

task SnvstomnvsAnnotateExomeWgsCalled {
    input {
        String pairName
        String chrom
        String finalChromVcfPath = "~{pairName}.mnv.final.v7.filtered.~{chrom}.vcf"
        File filteredOutFile
        Int memoryGb = 16
        Int diskSize = 4
    }

    command {
        python \
        /SNVsToMNVs_AnnotateExomeWgsCalled.py \
        -i ~{filteredOutFile} \
        -o ~{finalChromVcfPath}
    }

    output {
        File finalChromVcf = "~{finalChromVcfPath}"
    }

    runtime {
        mem: memoryGb + "G"
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/somatic_dna_tools@sha256:f281e73ddf515f3a5db6e766ef12fb2331dafa2b27c88e426764dcd8b4a7e19b"
        runtime_minutes: "90"
    }
}


task CompressVcf {
    input {
        File vcf
        String vcfCompressedPath = sub(basename(vcf), "$", ".gz")
        Int memoryGb = 4
        Int diskSize = (ceil( size(vcf, "GB") )  * 2 ) + 4
    }

    command {
        bgzip \
        -c \
        ~{vcf} \
        > ~{vcfCompressedPath} \
    }

    output {
        File vcfCompressed = "~{vcfCompressedPath}"
    }

    runtime {
        mem: memoryGb + "G"
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/samtools@sha256:963b0b2f24908832efab8ddccb7a7f3ba5dca9803bc099be7cf3a455766610fd"
        runtime_minutes: "500"
    }
}

task IndexVcf {
    input {
        File vcfCompressed
        String vcfCompressedPath = basename(vcfCompressed)
        String vcfIndexedPath = sub(basename(vcfCompressed), "$", ".tbi")
        Int threads = 4
        Int memoryGb = 8
        Int diskSize = (ceil( size(vcfCompressed, "GB") )  * 2 ) + 4
    }

    Int jvmHeap = memoryGb * 750  # Heap size in Megabytes. mem is in GB. (75% of mem)
    command {
        set -e -o pipefail

        cp \
        ~{vcfCompressed} \
        ~{vcfCompressedPath}

        gatk \
        IndexFeatureFile \
        --java-options "-Xmx~{jvmHeap}m -XX:ParallelGCThreads=4" \
        --input ~{vcfCompressedPath} \
        -O ~{vcfIndexedPath}
    }

    output {
        IndexedVcf vcfCompressedIndexed = object {
                vcf : "~{vcfCompressedPath}",
                index : "~{vcfIndexedPath}"
            }
    }

    runtime {
        mem: memoryGb + "G"
        cpus: threads
        cpu : threads
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/broadinstitute/gatk4@sha256:b3bde7bc74ab00ddce342bd511a9797007aaf3d22b9cfd7b52f416c893c3774c"
        runtime_minutes: "60"
    }
}


# Prep tasks

task PrepLancetConfirm {
    input {
        String pairName
        File callerVcf
        IndexedReference referenceFa
        # intermediate files
        String renameMetaVcfPath = sub(basename(callerVcf), "$", ".rename_metadata.vcf")
        String prepCallerVcfPath = sub(basename(callerVcf), "$", ".merge_prep.vcf")
        String renameVcfPath = sub(basename(callerVcf), "$", ".rename.vcf")
        String splitVcfPath = sub(basename(callerVcf), "$", ".split.vcf")
        String mnvVcfPath = sub(basename(callerVcf), "$", ".mnv.vcf")
        String tool
        String normalId
        String tumorId
        Int threads = 4
        Int memoryGb = 16
        Int diskSize = (ceil( size(callerVcf, "GB") )  * 2 ) + 4
    }

    command {
    set -e -o pipefail
    
    # task RenameMetadata
        python3 \
        /rename_metadata.py \
        ~{callerVcf} \
        ~{renameMetaVcfPath} \
        ~{tool}
    # task MergePrepSupport
        python3 \
        /merge_prep.py \
        --vcf ~{renameMetaVcfPath} \
        --out ~{prepCallerVcfPath} \
        --tool ~{tool} \
        --support
    # task RenameVcf
        python3 \
        /rename_vcf.py \
        ~{prepCallerVcfPath} \
        ~{renameVcfPath} \
        ~{normalId} \
        ~{tumorId} \
        ~{tool}
    # task confirmedIndexCompressVcf
        bgzip -c \
        ~{renameVcfPath} \
        > ~{renameVcfPath}.gz
        tabix \
        -p vcf \
        ~{renameVcfPath}.gz
    # task SplitMultiAllelic
        bcftools \
        norm \
        -m \
        -any \
        --threads ~{threads} \
        --no-version \
        -f ~{referenceFa.fasta} \
        -o ~{splitVcfPath} \
        ~{renameVcfPath}.gz
    # task SplitMnv {
        python3 \
        /split_mnv.py \
        ~{splitVcfPath} \
        ~{mnvVcfPath} \
        ~{tool}
    }

    output {
        File mnvVcf = "~{mnvVcfPath}"
    }

    runtime {
        mem: memoryGb + "G"
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-comp-s-fd4e/gcs_htslib_suite@sha256:47eba58683641905a58f31c05605c953dc1d888466e8d434a6ceff47b76df03e"
        runtime_minutes: "720"
    }
}


task RenameMetadata {
    input {
        String pairName
        File callerVcf
        String renameMetaVcfPath = sub(basename(callerVcf), "$", ".rename_metadata.vcf")
        String tool
        Int memoryGb = 16
        Int diskSize = (ceil( size(callerVcf, "GB") )  * 2 ) + 4
    }

    command {
        python \
        /rename_metadata.py \
        ~{callerVcf} \
        ~{renameMetaVcfPath} \
        ~{tool}
    }

    output {
        File renameMetaVcf = "~{renameMetaVcfPath}"
    }

    runtime {
        mem: memoryGb + "G"
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/somatic_dna_tools@sha256:f281e73ddf515f3a5db6e766ef12fb2331dafa2b27c88e426764dcd8b4a7e19b"
        runtime_minutes: "60"
    }
}

task MergePrepSupport {
    input {
        String pairName
        File renameMetaVcf
        String prepCallerVcfPath =  sub(basename(renameMetaVcf), ".rename_metadata.vcf$", ".merge_prep.vcf")
        String tool
        Int memoryGb = 16
        Int diskSize = (ceil( size(renameMetaVcf, "GB") )  * 2 ) + 4
    }

    command {
        python \
        /merge_prep.py \
        --vcf ~{renameMetaVcf} \
        --out ~{prepCallerVcfPath} \
        --tool ~{tool} \
        --support
    }

    output {
        File prepCallerVcf = "~{prepCallerVcfPath}"
    }

    runtime {
        mem: memoryGb + "G"
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/somatic_dna_tools@sha256:f281e73ddf515f3a5db6e766ef12fb2331dafa2b27c88e426764dcd8b4a7e19b"
        runtime_minutes: "60"
    }
}

task MergePrep {
    input {
        String pairName
        File renameMetaVcf
        String prepCallerVcfPath = sub(basename(renameMetaVcf, ".gz"), ".rename_metadata.vcf$", ".merge_prep.vcf")
        String tool
        Int memoryGb = 16
        Int diskSize = (ceil( size(renameMetaVcf, "GB") )  * 2 ) + 4
    }

    command {
        python \
        /merge_prep.py \
        --vcf ~{renameMetaVcf} \
        --out ~{prepCallerVcfPath} \
        --tool ~{tool}
    }

    output {
        File prepCallerVcf = "~{prepCallerVcfPath}"
    }

    runtime {
        mem: memoryGb + "G"
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/somatic_dna_tools@sha256:f281e73ddf515f3a5db6e766ef12fb2331dafa2b27c88e426764dcd8b4a7e19b"
        runtime_minutes: "60"
    }
}

task RenameVcf {
    input {
        File prepCallerVcf
        String pairName
        String renameVcfPath = sub(basename(prepCallerVcf, ".gz"), ".merge_prep.vcf$", ".rename.vcf")
        String normalId
        String tumorId
        String tool
        Int memoryGb = 16
        Int diskSize = (ceil( size(prepCallerVcf, "GB") )  * 2 ) + 4
    }

    command {
        python \
        /rename_vcf.py \
        ~{prepCallerVcf} \
        ~{renameVcfPath} \
        ~{normalId} \
        ~{tumorId} \
        ~{tool}
    }

    output {
        File renameVcf = "~{renameVcfPath}"
    }

    runtime {
        mem: memoryGb + "G"
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/somatic_dna_tools@sha256:f281e73ddf515f3a5db6e766ef12fb2331dafa2b27c88e426764dcd8b4a7e19b"
        runtime_minutes: "60"
    }
}

task RenameVcfPon {
    input {
        File prepCallerVcf
        String renameVcfPath = sub(basename(prepCallerVcf, ".gz"), ".vcf$", ".rename.vcf")
        String tumor
        String tool
        Int memoryGb = 16
        Int diskSize = (ceil( size(prepCallerVcf, "GB") )  * 2 ) + 4
    }

    command {
        python \
        /rename_vcf_pon.py \
        ~{prepCallerVcf} \
        ~{renameVcfPath} \
        ~{tumor} \
        ~{tool}
    }

    output {
        File renameVcf = "~{renameVcfPath}"
    }

    runtime {
        mem: memoryGb + "G"
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/somatic_dna_tools@sha256:f281e73ddf515f3a5db6e766ef12fb2331dafa2b27c88e426764dcd8b4a7e19b"
        runtime_minutes: "60"
    }
}

task SplitMultiAllelic {
    input {
        String pairName
        String splitVcfPath
        IndexedReference referenceFa
        IndexedVcf vcfCompressedIndexed
        Int threads = 8
        Int memoryGb = 16
        Int diskSize = (ceil( size(vcfCompressedIndexed.vcf, "GB") )  * 2 ) + 10
    }

    command {
        bcftools \
        norm \
        -m \
        -any \
        --threads ~{threads} \
        --no-version \
        -f ~{referenceFa.fasta} \
        -o ~{splitVcfPath} \
        ~{vcfCompressedIndexed.vcf}
    }

    output {
        File splitVcf = "~{splitVcfPath}"
    }

    runtime {
        mem: memoryGb + "G"
        cpus: threads
        cpu : threads
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/bcftools@sha256:d9b84254b8cc29fcae76b728a5a9a9a0ef0662ee52893d8d82446142876fb400"
        runtime_minutes: "90"
    }
}

task SplitMultiAllelicRegions {
    input {
        String pairName
        String splitVcfPath
        IndexedReference referenceFa
        IndexedVcf vcfCompressedIndexed
        Array[String]+ listOfChroms
        Int threads = 8
        Int memoryGb = 16
        Int diskSize = (ceil( size(vcfCompressedIndexed.vcf, "GB") )  * 3 ) + 10
    }

    command {
        bcftools \
        norm \
        -m \
        -any \
        --threads ~{threads} \
        --regions ~{sep="," listOfChroms} \
        --no-version \
        -f ~{referenceFa.fasta} \
        -o ~{splitVcfPath} \
        ~{vcfCompressedIndexed.vcf}
    }

    output {
        File sortedVcf = "~{splitVcfPath}"
    }

    runtime {
        mem: memoryGb + "G"
        cpus: threads
        cpu : threads
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/bcftools@sha256:d9b84254b8cc29fcae76b728a5a9a9a0ef0662ee52893d8d82446142876fb400"
        runtime_minutes: "90"
    }
}

task SplitMnv {
    input {
        File splitVcf
        String mnvVcfPath = sub(basename(splitVcf), ".split.vcf", ".split_mnvs.vcf")
        String tool
        Int memoryGb = 16
        Int diskSize = (ceil( size(splitVcf, "GB") )  * 2 ) + 4
    }

    command {
        python \
        /split_mnv.py \
        ~{splitVcf} \
        ~{mnvVcfPath} \
        ~{tool}
    }

    output {
        File mnvVcf = "~{mnvVcfPath}"
    }

    runtime {
        mem: memoryGb + "G"
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/somatic_dna_tools@sha256:f281e73ddf515f3a5db6e766ef12fb2331dafa2b27c88e426764dcd8b4a7e19b"
        runtime_minutes: "60"
    }
}

task RemoveContig {
    input {
        String mnvVcfPath
        String removeChromVcfPath = "~{mnvVcfPath}"
        File removeChromVcf
        String dollarSign = "$"
        Int memoryGb = 16
        Int diskSize = (ceil( size(removeChromVcf, "GB") )  * 2 ) + 4
    }

    command <<<
        set -e -o pipefail

        filename=$(basename -- "~{dollarSign}~{removeChromVcf}")
        extension="~{dollarSign}{filename##*.}"
        filename="~{dollarSign}{filename%.*}"

        if [[ ~{dollarSign}extension == gz ]]; then
            input_path=~{dollarSign}filename

            gunzip -c \
            ~{removeChromVcf} \
            > ~{dollarSign}input_path
        else
            input_path=~{removeChromVcf}
        fi

        python \
        /remove_contig.py \
        ~{dollarSign}input_path \
        ~{removeChromVcfPath}
    >>>

    output {
        File removeContigVcf = "~{removeChromVcfPath}"
    }

    runtime {
        mem: memoryGb + "G"
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/somatic_dna_tools@sha256:f281e73ddf515f3a5db6e766ef12fb2331dafa2b27c88e426764dcd8b4a7e19b"
        runtime_minutes: "60"
    }
}

task Gatk4MergeSortVcf {
    input {
        String sortedVcfPath
        Array[File] tempVcfs
        IndexedReference referenceFa
        Boolean gzipped = false
        String suffix = if gzipped then ".tbi" else ".idx"
        Int threads = 4
        Int memoryGb = 8
        Int diskSize = 10
    }

    Int jvmHeap = memoryGb * 750  # Heap size in Megabytes. mem is in GB. (75% of mem)

    command {
        gatk \
        SortVcf \
        --java-options "-Xmx~{jvmHeap}m -XX:ParallelGCThreads=4" \
        -SD ~{referenceFa.dict} \
        -I ~{sep=" -I " tempVcfs} \
        -O ~{sortedVcfPath}
    }

    output {
        IndexedVcf sortedVcf = object {
            vcf : "~{sortedVcfPath}",
            index : "~{sortedVcfPath}~{suffix}"
        }
    }

    runtime {
        mem: memoryGb + "G"
        cpus: threads
        cpu : threads
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/broadinstitute/gatk4@sha256:b3bde7bc74ab00ddce342bd511a9797007aaf3d22b9cfd7b52f416c893c3774c"
        runtime_minutes: "500"
    }
}

# Merge Callers Section
task MergeCenterCallers {
    input {
        String pairName
        # input files
        Array[IndexedVcf] allVcfCompressed
        Array[File] allVcfCompressedList
        String mergedVcfPath = "~{pairName}.merged.v7.vcf"
        String infoTags = "centers_called_by:join,num_centers:sum,called_by:join,MNV_ID:join"
        Int threads = 4
        Int memoryGb = 16
        Int diskSize = 20
    }
    command {
    set -e -o pipefail
    # task MergeCenterCallers
        bcftools \
        merge \
        --force-samples \
        --no-version \
        --threads ~{threads} \
        -F x \
        -m none \
        -o ~{mergedVcfPath} \
        -i ~{infoTags} \
        ~{sep=" " allVcfCompressedList}
    }
    output {
        File mergedVcf = "~{mergedVcfPath}"
    }

    runtime {
        mem: memoryGb + "G"
        disks: "local-disk " + diskSize + " HDD"
        cpus: threads
        cpu : threads
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-comp-s-fd4e/gcs_htslib_suite@sha256:47eba58683641905a58f31c05605c953dc1d888466e8d434a6ceff47b76df03e"
        runtime_minutes: "720"
    }
}

task CountCenterCallers  {
    input {
        String pairName
        String chrom
        String tumorId
        String normalId
        Bam normalFinalBam
        Bam tumorFinalBam
        File mergedVcf
        String mergedVcfCompressedPath = "~{pairName}.merged_supported.v7.vcf.gz"
        String mergedVcfIndexedPath = "~{pairName}.merged_supported.v7.vcf.gz.tbi"
        # intermediate files
        String columnChromVcfPath = "~{pairName}.single_column.v7.~{chrom}.vcf"
        # input files
        Int threads = 1
        Int memoryGb = 16
        Int diskSize = 20
        String preCountsChromVcfPath = "~{pairName}.pre_count.v7.~{chrom}.vcf"
        String countsChromVcfPath = "~{pairName}.final.v7.~{chrom}.vcf"
        String finalChromVcfPath = "~{chrom}.mnv.final.v7.filtered.vcf"
        File mergeCenterColumnsScript = "gs://nygc-comp-s-fd4e-input/somatic/hcmi/tools/merge_columns.py"
        File snvsToMnvs = "gs://nygc-comp-s-fd4e-input/somatic/hcmi/tools/SNVsToMNVs_CountsBasedFilter_AnnotateHighConfCenter.py"
        File Classes = "gs://nygc-comp-s-fd4e-input/somatic/hcmi/tools/Classes.py"
        File AddNygcAlleleCountsToVcf = "gs://nygc-comp-s-fd4e-input/somatic/hcmi/tools/add_nygc_allele_counts_to_vcf_center.py"
#         File mergeCenterColumnsScript = "pipelines/snv-indel-concordance/somatic_dna_tools/merge_columns.py"
#         File snvsToMnvs = "pipelines/snv-indel-concordance/somatic_dna_tools/SNVsToMNVs_CountsBasedFilter_AnnotateHighConfCenter.py"
#         File Classes = "pipelines/snv-indel-concordance/somatic_dna_tools/Classes.py"
#         File AddNygcAlleleCountsToVcf = "pipelines/snv-indel-concordance/somatic_dna_tools/add_nygc_allele_counts_to_vcf_center.py"
    }

    command {
        set -e -o pipefail
    # MergeCenterColumns
        python ~{mergeCenterColumnsScript} \
        ~{mergedVcf} \
        ~{columnChromVcfPath} \
        ~{tumorId} \
        ~{normalId}
    # AddNygcAlleleCountsToVcf
        python3 \
        ~{AddNygcAlleleCountsToVcf} \
        -t ~{tumorFinalBam.bam} \
        -n ~{normalFinalBam.bam} \
        -v ~{columnChromVcfPath} \
        -b 10 \
        -m 10 \
        -o ~{countsChromVcfPath}
    # AddFinalAlleleCountsToVcf (skipped)
    # SnvstomnvsCountsbasedfilterAnnotatehighconfCenter
        python \
        ~{snvsToMnvs} \
        -i ~{countsChromVcfPath} \
        -o ~{finalChromVcfPath}
    }

    output {
        File mergedChromVcf = "~{columnChromVcfPath}"
        File finalChromVcf = "~{finalChromVcfPath}"
    }

    runtime {
        mem: memoryGb + "G"
        disks: "local-disk " + diskSize + " HDD"
        cpus: threads
        cpu : threads
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-comp-s-fd4e/gcs_htslib_suite@sha256:47eba58683641905a58f31c05605c953dc1d888466e8d434a6ceff47b76df03e"
        runtime_minutes: "720"
    }
}

task MergeCallersGetCandidate {
    input {
        String pairName
        String mergedVcfPath = "~{pairName}.merged_supported.v7.vcf"
        String mergedVcfCompressedPath = "~{pairName}.merged_supported.v7.vcf.gz"
        String mergedVcfIndexedPath = "~{pairName}.merged_supported.v7.vcf.gz.tbi"
        String candidateBedPrefixPath = "~{pairName}.candidate.split."
        # intermediate files
        String startVcfPath = "~{pairName}.start.merged.v7.vcf"
        String candidateVcfPath = "~{pairName}.candidate.merged.v7.vcf"
        String candidateBedPath = "~{pairName}.candidate.bed"
        Array[String] candidateBedPaths
        # input files
        File intervalListBed
        Array[IndexedVcf] allVcfCompressed
        Array[File] allVcfCompressedList
        Int minSplits
        Int threads = 8
        Int memoryGb = 16
        Int diskSize = 20
        String additionalSuffix = ".bed"
    }

    command {
    set -e -o pipefail
    # task MergeCallers
        bcftools \
        merge \
        --force-samples \
        --no-version \
        --threads ~{threads} \
        -f PASS,SUPPORT \
        -F x \
        -m none \
        -o ~{mergedVcfPath} \
        -i called_by:join,num_callers:sum,MNV_ID:join,supported_by:join \
        ~{sep=" " allVcfCompressedList}
    # task StartCandidates
        bedtools \
        intersect \
        -header \
        -a ~{mergedVcfPath} \
        -b ~{intervalListBed} \
        -v \
        > ~{startVcfPath}
    # task GetCandidates
        python3 \
        /get_candidates.py \
        ~{startVcfPath} \
        ~{candidateVcfPath}
    # task VcfToBed
        python3 \
        /vcf_to_bed.py \
        ~{candidateVcfPath} \
        | bedtools \
        merge \
        > ~{candidateBedPath}
        
        split \
        --number=l/~{minSplits} \
        -d \
        -a2 \
        --additional-suffix=~{additionalSuffix} \
        ~{candidateBedPath} \
        ~{candidateBedPrefixPath}
        
        bgzip -c \
        ~{mergedVcfPath} \
        > ~{mergedVcfCompressedPath}
        tabix \
        -p vcf \
        ~{mergedVcfCompressedPath}
        
        bgzip -c \
        ~{candidateVcfPath} \
        > ~{candidateVcfPath}.gz
        tabix \
        -p vcf \
        ~{candidateVcfPath}.gz
    }

    output {
        Array[File?] candidateBeds = candidateBedPaths
        File mergedVcf = "~{mergedVcfPath}"
        IndexedVcf mergedVcfCompressedIndexed = object {
            vcf : "~{mergedVcfCompressedPath}",
            index : "~{mergedVcfIndexedPath}"
        }
        IndexedVcf candidateVcfCompressedIndexed = object {
            vcf : "~{candidateVcfPath}.gz",
            index : "~{candidateVcfPath}.gz.tbi"
        }
    }

    runtime {
        mem: memoryGb + "G"
        disks: "local-disk " + diskSize + " HDD"
        cpus: threads
        cpu : threads
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-comp-s-fd4e/gcs_htslib_suite@sha256:47eba58683641905a58f31c05605c953dc1d888466e8d434a6ceff47b76df03e"
        runtime_minutes: "720"
    }
}

task MergeCallersAll {
    input {
        String pairName
        String mergedChromVcfPath = "~{pairName}.merged_supported.v7.vcf"
        Array[IndexedVcf] allVcfCompressed
        Array[File] allVcfCompressedList
        Int threads = 8
        Int memoryGb = 16
        Int diskSize = 20
    }

    command {
        bcftools \
        merge \
        --force-samples \
        --no-version \
        --threads ~{threads} \
        -f PASS,SUPPORT \
        -F x \
        -m none \
        -o ~{mergedChromVcfPath} \
        -i called_by:join,num_callers:sum,MNV_ID:join,supported_by:join \
        ~{sep=" " allVcfCompressedList}
    }

    output {
        File mergedChromVcf = "~{mergedChromVcfPath}"
    }

    runtime {
        mem: memoryGb + "G"
        cpus: threads
        cpu : threads
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/bcftools@sha256:d9b84254b8cc29fcae76b728a5a9a9a0ef0662ee52893d8d82446142876fb400"
        runtime_minutes: "60"
    }
}

task MergeCallers {
    input {
        String chrom
        String pairName
        String mergedChromVcfPath = "~{pairName}.merged_supported.v7.~{chrom}.vcf"
        Array[IndexedVcf] allVcfCompressed
        Array[File] allVcfCompressedList
        Int threads = 8
        Int memoryGb = 16
        Int diskSize = 20
    }

    command {
        bcftools \
        merge \
        -r ~{chrom} \
        --force-samples \
        --no-version \
        --threads ~{threads} \
        -f PASS,SUPPORT \
        -F x \
        -m none \
        -o ~{mergedChromVcfPath} \
        -i called_by:join,num_callers:sum,MNV_ID:join,supported_by:join \
        ~{sep=" " allVcfCompressedList}
    }

    output {
        File mergedChromVcf = "~{mergedChromVcfPath}"
    }

    runtime {
        mem: memoryGb + "G"
        cpus: threads
        cpu : threads
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/bcftools@sha256:d9b84254b8cc29fcae76b728a5a9a9a0ef0662ee52893d8d82446142876fb400"
        runtime_minutes: "60"
    }
}


task StartCandidates {
    input {
        String pairName
        String chrom
        String startChromVcfPath = "~{pairName}.start.merged.v7.~{chrom}.vcf"
        File intervalListBed
        File mergedChromVcf
        Int memoryGb = 16
        Int diskSize = 20
    }

    command {
        bedtools \
        intersect \
        -header \
        -a ~{mergedChromVcf} \
        -b ~{intervalListBed} \
        -v \
        > ~{startChromVcfPath}
    }

    output {
        File startChromVcf = "~{startChromVcfPath}"
    }

    runtime {
        mem: memoryGb + "G"
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/bedtools@sha256:9e737f5c96c00cf3b813d419d7a7b474c4013c9aa9dfe704eb36417570c6474e"
        runtime_minutes: "60"
    }
}

task GetCandidates {
    input {
        String pairName
        String chrom
        String candidateChromVcfPath = "~{pairName}.candidate.merged.v7.~{chrom}.vcf"
        File startChromVcf
        Int memoryGb = 16
        Int diskSize = (ceil( size(startChromVcf, "GB") )  * 2 ) + 4
    }

    command {
        python \
        /get_candidates.py \
        ~{startChromVcf} \
        ~{candidateChromVcfPath}
    }

    output {
        File candidateChromVcf = "~{candidateChromVcfPath}"
    }

    runtime {
        mem: memoryGb + "G"
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/somatic_dna_tools@sha256:f281e73ddf515f3a5db6e766ef12fb2331dafa2b27c88e426764dcd8b4a7e19b"
        runtime_minutes: "60"
    }
}

task VcfToBed {
    input {
        String pairName
        String chrom
        String candidateChromBedPath = "~{pairName}.candidate.merged.v7.~{chrom}.bed"
        File candidateChromVcf
        Int memoryGb = 16
        Int diskSize = (ceil( size(candidateChromVcf, "GB") )  * 2 ) + 4
    }

    command {
        set -e -o pipefail

        python \
        /vcf_to_bed.py \
        ~{candidateChromVcf} \
        | bedtools \
        merge \
        > ~{candidateChromBedPath}
    }

    output {
        File candidateChromBed = "~{candidateChromBedPath}"
    }

    runtime {
        mem: memoryGb + "G"
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/somatic_dna_tools@sha256:f281e73ddf515f3a5db6e766ef12fb2331dafa2b27c88e426764dcd8b4a7e19b"
        runtime_minutes: "90"
    }
}

task LancetConfirm {
    input {
        String pairName
        String chrom
        String lancetChromVcfPath = "~{pairName}.lancet.merged.v7.~{chrom}.vcf"
        IndexedReference referenceFa
        Bam normalFinalBam
        File candidateChromBed
        Bam tumorFinalBam
        Int threads = 8
        Int memoryGb = 16
        Int diskSize = (ceil( size(tumorFinalBam.bam, "GB") + size(normalFinalBam.bam, "GB")) ) + 10
    }

    command {
        set -e -o pipefail

        mkdir ~{chrom}

        cd ~{chrom}

        lancet \
        --normal ~{normalFinalBam.bam} \
        --tumor ~{tumorFinalBam.bam} \
        --bed ~{candidateChromBed} \
        --ref ~{referenceFa.fasta} \
        --min-k 11 \
        --low-cov 1 \
        --min-phred-fisher 5 \
        --min-strand-bias 1 \
        --min-alt-count-tumor 3 \
        --min-vaf-tumor 0.04 \
        --padding 250 \
        --window-size 2000 \
        --num-threads ~{threads} \
        > ~{lancetChromVcfPath}
    }

    output {
        File lancetChromVcf = "~{chrom}/~{lancetChromVcfPath}"
    }

    runtime {
        mem: memoryGb + "G"
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/lancet@sha256:25169d34b41de9564e03f02ebcbfb4655cf536449592b0bd58773195f9376e61"
        runtime_minutes: "1000"
    }
}

task IntersectVcfs {
    input {
        Int threads = 8
        Int memoryGb = 16
        Int diskSize = 4
        String pairName
        String chrom
        String vcfConfirmedCandidatePath = "~{pairName}.confirmed_lancet.merged.v7.~{chrom}.vcf"
        IndexedVcf vcfCompressedLancet
        IndexedVcf vcfCompressedCandidate
    }

    command {
        bcftools \
        isec \
        -w 1 \
        -c none \
        -n =2 \
        --threads ~{threads} \
        ~{vcfCompressedLancet.vcf} \
        ~{vcfCompressedCandidate.vcf} \
        > ~{vcfConfirmedCandidatePath}
    }

    output {
        File vcfConfirmedCandidate = "~{vcfConfirmedCandidatePath}"
    }

    runtime {
        mem: memoryGb + "G"
        cpus: threads
        cpu : threads
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/bcftools@sha256:d9b84254b8cc29fcae76b728a5a9a9a0ef0662ee52893d8d82446142876fb400"
        runtime_minutes: "60"
    }
}

task PostProcessMerged {
    input {
        String pairName
        String chrom
        File ponFile
        # intermediate files
        String columnChromVcfPath = "~{pairName}.single_column.v7.~{chrom}.vcf"
        String columnChromBedPath = "~{pairName}.single_column.v7.~{chrom}.bed"
        String preCountsChromVcfPath = "~{pairName}.pre_count.v7.~{chrom}.vcf"
        String countsChromVcfPath = "~{pairName}.final.v7.~{chrom}.vcf"
        String ponOutFilePath = "~{pairName}.final.v7.pon.~{chrom}.vcf"
        String finalChromVcfPath = "~{pairName}.mnv.final.v7.filtered.~{chrom}.vcf"
        String tumorId
        String normalId
        String tumorBamSlicePath = "~{tumorId}.sliced.bam"
        String tumorBaiSlicePath = "~{tumorId}.sliced.bai"
        String normalBamSlicePath = "~{normalId}.sliced.bam"
        String normalBaiSlicePath = "~{normalId}.sliced.bai"
        File supportedChromVcf
        Bam normalFinalBam
        Bam tumorFinalBam
        Int memoryGb = 60
        Int threads = 4
        Int diskSize = 20
        
        # needed for work on cloud without localizing
        File? serviceAccountKey
        String? gcpProject
    }

    command {
    set -e -o pipefail
    
    serviceAccountKey=~{serviceAccountKey}
    if [ -f "$serviceAccountKey" ]; then
        export GOOGLE_APPLICATION_CREDENTIALS=~{serviceAccountKey}
        export GOOGLE_CLOUD_PROJECT=${gcpProject}
        # expires in 60 min
        export GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token`
        export GCS_REQUESTER_PAYS_PROJECT=~{gcpProject}
    fi
        
    # MergeColumns
        python3 \
        /merge_columns.py \
        ~{supportedChromVcf} \
        ~{columnChromVcfPath} \
        ~{tumorId} \
        ~{normalId}
        
    # task VcfToBed
        python3 \
        /vcf_to_bed.py \
        ~{columnChromVcfPath} \
        | bedtools \
        merge \
        > ~{columnChromBedPath}
    # SliceBam Tumor
        samtools view \
        --threads ~{threads} \
        -h \
        -b \
        -L ~{columnChromBedPath} \
        ~{tumorFinalBam.bam} \
        > ~{tumorBamSlicePath}
        
        samtools \
        index \
        -@ ~{threads} \
        ~{tumorBamSlicePath} \
        ~{tumorBaiSlicePath}

    # SliceBam Normal
        samtools view \
        --threads ~{threads} \
        -h \
        -b \
        -L ~{columnChromBedPath} \
        ~{normalFinalBam.bam} \
        > ~{normalBamSlicePath}
        
        samtools \
        index \
        -@ ~{threads} \
        ~{normalBamSlicePath} \
        ~{normalBaiSlicePath}
    
    # AddNygcAlleleCountsToVcf
        python3 \
        /add_nygc_allele_counts_to_vcf.py \
        -t ~{tumorBamSlicePath} \
        -n ~{normalBamSlicePath} \
        -v ~{columnChromVcfPath} \
        -b 10 \
        -m 10 \
        -o ~{preCountsChromVcfPath}
    # AddFinalAlleleCountsToVcf
        python3 \
        /add_final_allele_counts_to_vcf.py \
        -v ~{preCountsChromVcfPath} \
        -o ~{countsChromVcfPath}
    # FilterPon
        python3 \
        /filter_pon.py \
        --bed ~{ponFile} \
        --vcf ~{countsChromVcfPath} \
        --out ~{ponOutFilePath}
    # SnvstomnvsCountsbasedfilterAnnotatehighconf
        python3 \
        /SNVsToMNVs_CountsBasedFilter_AnnotateHighConf.py \
        -i ~{ponOutFilePath} \
        -o ~{finalChromVcfPath}
    }
    
    output {
        File finalChromVcf = "~{finalChromVcfPath}"
    }

    runtime {
        mem: memoryGb + "G"
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-comp-s-fd4e/gcs_htslib_suite@sha256:47eba58683641905a58f31c05605c953dc1d888466e8d434a6ceff47b76df03e"
        runtime_minutes: "1440"
    }
    parameter_meta {
        normalFinalBam: {
            description: "Input BAM filename",
            category: "required",
            localization_optional: true
        }
        tumorFinalBam: {
            description: "Input BAM filename",
            category: "required",
            localization_optional: true
        }
    }
}
    
task MergeColumns {
    input {
        String pairName
        String chrom
        String columnChromVcfPath = "~{pairName}.single_column.v7.~{chrom}.vcf"
        String tumorId
        String normalId
        File supportedChromVcf
        Int memoryGb = 16
        Int diskSize = 4
    }

    command {
        python \
        /merge_columns.py \
        ~{supportedChromVcf} \
        ~{columnChromVcfPath} \
        ~{tumorId} \
        ~{normalId}
    }

    output {
        File columnChromVcf = "~{columnChromVcfPath}"
    }

    runtime {
        mem: memoryGb + "G"
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/somatic_dna_tools@sha256:f281e73ddf515f3a5db6e766ef12fb2331dafa2b27c88e426764dcd8b4a7e19b"
        runtime_minutes: "60"
    }
}

task MergeColumnsPon {
    input {
        String chrom
        String tumor
        String columnChromVcfPath = "~{tumor}.single_column.v7.~{chrom}.vcf"
        File supportedChromVcf
        Int memoryGb = 16
        Int diskSize = 4
    }

    command {
        python \
        /merge_columns_pon.py \
        ~{supportedChromVcf} \
        ~{columnChromVcfPath} \
        ~{tumor}
    }

    output {
        File columnChromVcf = "~{columnChromVcfPath}"
    }

    runtime {
        mem: memoryGb + "G"
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/somatic_dna_tools@sha256:f281e73ddf515f3a5db6e766ef12fb2331dafa2b27c88e426764dcd8b4a7e19b"
        runtime_minutes: "60"
    }
}

task AddNygcAlleleCountsToVcf {
    input {
        String pairName
        String chrom
        String preCountsChromVcfPath = "~{pairName}.pre_count.v7.~{chrom}.vcf"
        Bam normalFinalBam
        Bam tumorFinalBam
        File columnChromVcf
        Int memoryGb = 40
        Int diskSize = ceil( size(tumorFinalBam.bam, "GB") + size(normalFinalBam.bam, "GB")) + 20
    }

    command {
        python \
        /add_nygc_allele_counts_to_vcf.py \
        -t ~{tumorFinalBam.bam} \
        -n ~{normalFinalBam.bam} \
        -v ~{columnChromVcf} \
        -b 10 \
        -m 10 \
        -o ~{preCountsChromVcfPath}
    }

    output {
        File preCountsChromVcf = "~{preCountsChromVcfPath}"
    }

    runtime {
        mem: memoryGb + "G"
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/somatic_dna_tools@sha256:f281e73ddf515f3a5db6e766ef12fb2331dafa2b27c88e426764dcd8b4a7e19b"
        runtime_minutes: "500"
    }
}

task AddFinalAlleleCountsToVcf {
    input {
        String pairName
        String chrom
        String countsChromVcfPath = "~{pairName}.final.v7.~{chrom}.vcf"
        File preCountsChromVcf
        Int memoryGb = 16
        Int diskSize = 4
    }

    command {
        python \
        /add_final_allele_counts_to_vcf.py \
        -v ~{preCountsChromVcf} \
        -o ~{countsChromVcfPath} \
    }

    output {
        File countsChromVcf = "~{countsChromVcfPath}"
    }

    runtime {
        mem: memoryGb + "G"
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/somatic_dna_tools@sha256:f281e73ddf515f3a5db6e766ef12fb2331dafa2b27c88e426764dcd8b4a7e19b"
        runtime_minutes: "60"
    }
}

task FilterPon {
    input {
        String chrom
        String pairName
        String ponOutFilePath = "~{pairName}.pon.final.v7.~{chrom}.vcf"
        File countsChromVcf
        File ponFile
        Int memoryGb = 60
        Int diskSize = 16
    }

    command {
        python \
        /filter_pon.py \
        --bed ~{ponFile} \
        --chrom ~{chrom} \
        --vcf ~{countsChromVcf} \
        --out ~{ponOutFilePath} \
    }

    output {
        File ponOutFile = "~{ponOutFilePath}"
    }

    runtime {
        mem: memoryGb + "G"
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/somatic_dna_tools@sha256:f281e73ddf515f3a5db6e766ef12fb2331dafa2b27c88e426764dcd8b4a7e19b"
        runtime_minutes: "90"
    }
}

task SnvstomnvsCountsbasedfilterAnnotatehighconf {
    input {
        String pairName
        String chrom
        String finalChromVcfPath = "~{pairName}.mnv.final.v7.filtered.~{chrom}.vcf"
        File filteredOutFile
        Int memoryGb = 16
        Int diskSize = 4
    }

    command {
        python \
        /SNVsToMNVs_CountsBasedFilter_AnnotateHighConf.py \
        -i ~{filteredOutFile} \
        -o ~{finalChromVcfPath} \
    }

    output {
        File finalChromVcf = "~{finalChromVcfPath}"
    }

    runtime {
        mem: memoryGb + "G"
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/somatic_dna_tools@sha256:f281e73ddf515f3a5db6e766ef12fb2331dafa2b27c88e426764dcd8b4a7e19b"
        runtime_minutes: "60"
    }
}
