version 1.0

import "../wdl_structs.wdl"

task FineGrainedCov {
    input {
        Bam finalBam
        Int windowSize = 1000
        String sampleId
        String analysisId
        String distMosdepthRegionAutoPath = "~{sampleId}.mosdepth.region.dist.txt"
        String summaryMosdepthAutoPath = "~{sampleId}.mosdepth.summary.txt"
        String regionsMosdepthAutoPath = "~{sampleId}.regions.bed.gz"
        String regionsMosdepthIndexAutoPath = "~{sampleId}.regions.bed.gz.csi"
        
        String distMosdepthRegionPath = "~{analysisId}.mosdepth.region.dist.txt"
        String summaryMosdepthPath = "~{analysisId}.mosdepth.summary.txt"
        String regionsMosdepthPath = "~{analysisId}.regions.bed.gz"
        String regionsMosdepthIndexPath = "~{analysisId}.regions.bed.gz.csi"
        # resources
        Int diskSize
        Int threads = 4
        Int memoryGb = 16
    }

    command {
        mosdepth \
        --by ~{windowSize} \
        --no-per-base \
        --threads ~{threads} \
        ~{sampleId} \
        ~{finalBam.bam}
        cp ~{distMosdepthRegionAutoPath} \
        ~{distMosdepthRegionAutoPath}.cp
        cp ~{distMosdepthRegionAutoPath}.cp \
        ~{distMosdepthRegionPath}
        
        cp ~{summaryMosdepthAutoPath} \
        ~{summaryMosdepthAutoPath}.cp
        cp ~{summaryMosdepthAutoPath}.cp \
        ~{summaryMosdepthPath}
        
        cp ~{regionsMosdepthAutoPath} \
        ~{regionsMosdepthAutoPath}.cp
        cp ~{regionsMosdepthAutoPath}.cp \
        ~{regionsMosdepthPath}
        
        cp ~{regionsMosdepthIndexAutoPath} \
        ~{regionsMosdepthIndexAutoPath}.cp
        cp ~{regionsMosdepthIndexAutoPath}.cp \
        ~{regionsMosdepthIndexPath}
    }

    output {
        File distMosdepthRegion = "~{distMosdepthRegionPath}"
        File summaryMosdepth = "~{summaryMosdepthPath}"
        File regionsMosdepth = "~{regionsMosdepthPath}"
        File regionsMosdepthIndex = "~{regionsMosdepthIndexPath}"
    }

    runtime {
        docker: "gcr.io/nygc-comp-s-fd4e/mosdepth@sha256:4e634382de98004957a2f50739492b2c8d57cd17584ff920b5c48036edd8dffe"
        disks: "local-disk " + diskSize + " HDD"
        memory: memoryGb + "GB"
        mem: memoryGb + "G"
        cpu: threads
        cpus: threads
        runtime_minutes: "100"
    }
}

task MakeVariantCram {
    input {
        Bam finalBam
        # mainVcf for SNVs, INDELs
        File mainVcf
        File cnvAnnotatedFinalBed
        File svFinalBedPe
        String sampleId
        String mainBedPath = "~{sampleId}.snv.indel.bed"
        String svChromABedPath = "~{sampleId}.sv.chromA.bed"
        String svChromBBedPath = "~{sampleId}.sv.chromB.bed"
        String cnvBedChromAPath = "~{sampleId}.cnv.chromA.bed"
        String cnvBedChromBPath = "~{sampleId}.cnv.chromB.bed"
        String featuresBedPath = "~{sampleId}.features.bed"
        String features500BedPath = "~{sampleId}.features.500.bed"
        String variantCramPath = "~{sampleId}.variantRegions.cram"
        String variantCraiPath = "~{sampleId}.variantRegions.crai"
        IndexedReference referenceFa
        Int slopSizeBp = 500
        String dollarSign = "$"
        # resources
        Int diskSize
        Int threads = 8
        Int memoryGb = 16
        
        # needed for work on cloud without localizing
        File? serviceAccountKey
        String? gcpProject
    }

    command <<<
        set -e -o pipefail
        
        cut -f1,2 \
        ~{referenceFa.index} \
        > sizes_genome.tsv
        
        python3 \
        /vcf_to_bed.py \
        ~{mainVcf} \
        | cut \
        -f 1,2,3 \
        | sort \
        -k 1,1 \
        -k2,2n \
        | bedtools \
        merge \
        > ~{mainBedPath}
        
        grep -v "^#" \
        ~{svFinalBedPe} \
        | cut \
        -f 1,2,3 \
        | sort \
        -k 1,1 \
        -k2,2n \
        | bedtools \
        merge \
        > ~{svChromABedPath}
        
        grep -v "^#" \
        ~{svFinalBedPe} \
        | cut \
        -f 4,5,6 \
        | sort \
        -k 1,1 \
        -k2,2n \
        | bedtools \
        merge \
        > ~{svChromBBedPath}
        
        grep -v "^#" \
        ~{cnvAnnotatedFinalBed} \
        | cut \
        -f 1,2 \
        | awk \
        'BEGIN{FS="\t"; OFS=FS}{print ~{dollarSign}1, ~{dollarSign}2, ~{dollarSign}2+1}' \
        | sort \
        -k 1,1 \
        -k2,2n \
        | bedtools \
        merge \
        > ~{cnvBedChromAPath}
        
        grep -v "^#" \
        ~{cnvAnnotatedFinalBed} \
        | cut \
        -f 1,3 \
        | awk \
        'BEGIN{FS="\t"; OFS=FS}{print ~{dollarSign}1, ~{dollarSign}2, ~{dollarSign}2+1}' \
        | sort \
        -k 1,1 \
        -k2,2n \
        | bedtools \
        merge \
        > ~{cnvBedChromBPath}
        
        cat \
        ~{mainBedPath} \
        ~{svChromABedPath} \
        ~{svChromBBedPath} \
        ~{cnvBedChromAPath} \
        ~{cnvBedChromBPath} \
        | sort \
        -k 1,1 \
        -k2,2n \
        | bedtools \
        merge \
        > ~{featuresBedPath}
        
        bedtools \
        slop \
        -i ~{featuresBedPath} \
        -g sizes_genome.tsv \
        -b ~{slopSizeBp} \
        > ~{features500BedPath}

        serviceAccountKey=~{serviceAccountKey}
        if [ -f "$serviceAccountKey" ]; then
            export GOOGLE_APPLICATION_CREDENTIALS=~{serviceAccountKey}
            export GOOGLE_CLOUD_PROJECT=${gcpProject}
            # expires in 60 min
            export GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token`
            export GCS_REQUESTER_PAYS_PROJECT=~{gcpProject}
        fi
        
        samtools view \
        --threads ~{threads} \
        -h \
        --cram \
        --reference ~{referenceFa.fasta} \
        -L ~{features500BedPath} \
        ~{finalBam.bam} \
        > ~{variantCramPath}
        
        samtools \
        index \
        -@ ~{threads} \
        ~{variantCramPath} \
        ~{variantCraiPath}
    >>>

    output {
        Cram variantCram = object {
                cram : "~{variantCramPath}",
                cramIndex : "~{variantCraiPath}"
             }
    }
    runtime {
        docker: "gcr.io/nygc-comp-s-fd4e/gcs_htslib_suite@sha256:47eba58683641905a58f31c05605c953dc1d888466e8d434a6ceff47b76df03e"
        disks: "local-disk " + diskSize + " HDD"
        memory: memoryGb + "GB"
        mem: memoryGb + "G"
        cpu: threads
        cpus: threads
        runtime_minutes: "600"
        maxRetries : 4
    }
    parameter_meta {
        finalBam: {
            description: "Input BAM filename",
            category: "required",
            localization_optional: true
        }
    }
}

task BamQcCheck {
    input {
        File wgsMetricsFile
        Float expectedCoverage
    }

    command {
        python /check_wgs_coverage.py -m ~{wgsMetricsFile} -e ~{expectedCoverage}
    }

    output {
        Boolean coveragePass = read_boolean(stdout())
    }

    runtime {
        docker: "gcr.io/nygc-public/workflow_utils@sha256:40fa18ac3f9d9f3b9f037ec091cb0c2c26ad6c7cb5c32fb16c1c0cf2a5c9caea"
        runtime_minutes: "60"
    }
}

task SomaticQcCheck {
    # Check coverage, contamination and concordance in one go.
    input {
        String pairId
        File tumorWgsMetricsFile
        File normalWgsMetricsFile
        Float tumorExpectedCoverage
        Float normalExpectedCoverage
        File concordanceFile
        File contaminationFile
        Float minConcordance = 95.0
        Float maxContamination = 0.99
        Boolean tumorSkipCoverageCheck = false
        Boolean normalSkipCoverageCheck = false
    }

    command {

        python /check_somatic_qc.py \
           --tumor_metrics_file ~{tumorWgsMetricsFile} \
           --tumor_expected_coverage ~{tumorExpectedCoverage} \
           --normal_metrics_file ~{normalWgsMetricsFile} \
           --normal_expected_coverage ~{normalExpectedCoverage} \
           --concordance_file ~{concordanceFile} \
           --contamination_file ~{contaminationFile} \
           --min_concordance ~{minConcordance} \
           --max_contamination ~{maxContamination} \
           ${if tumorSkipCoverageCheck then "--skip_tumor_coverage" else " "} \
           ${if normalSkipCoverageCheck then "--skip_normal_coverage" else " "}

        mv qc_summary.txt ~{pairId}.qc_summary.txt

    }

    output {
        Boolean qcPass = read_boolean(stdout())
        Boolean costPass = read_boolean(stdout())
        File qcCheckReport = "~{pairId}.qc_summary.txt"
    }

    runtime {
        docker: "gcr.io/nygc-public/workflow_utils@sha256:40fa18ac3f9d9f3b9f037ec091cb0c2c26ad6c7cb5c32fb16c1c0cf2a5c9caea"
        runtime_minutes: "60"
    }
}

task TestMutect2Rate {
    input {
        Int diskSize = 1
        Int memoryGb = 4
        String dollarSign = "$"
        Int mutect2RateThreshold = 1000
        File mutect2Rate
        String mutect2RatePassedPath = "mutect2.rate.passed"
    }

    command {
        number=~{mutect2RateThreshold}
        value=~{dollarSign}(cat ~{mutect2Rate})
        if [ "~{dollarSign}value" -lt "~{dollarSign}number" ]; then
            echo "The value in file.txt is less than $number" > ~{mutect2RatePassedPath}
        else
            echo "The value in ~{dollarSign}value is greater than or equal to ~{dollarSign}number"
        fi
    }

    output {
        File? mutect2RatePassed = "~{mutect2RatePassedPath}"
    }

    runtime {
        mem: memoryGb + "G"
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/somatic_dna_tools@sha256:f281e73ddf515f3a5db6e766ef12fb2331dafa2b27c88e426764dcd8b4a7e19b"
    }
}

task GetMutect2RateReadless {
    input {
        Int diskSize = 1
        Int memoryGb = 4
        Array[File] mutect2Logs
        String mutect2RateTablePath
        String mutect2RatePath = "mutect2.rate.txt"
    }

    command {
        python \
        /get_mutect_rate.py \
        ~{sep=" " mutect2Logs} \
        ~{mutect2RateTablePath} \
        > ~{mutect2RatePath}
    }

    output {
        File mutect2RateTable = "~{mutect2RateTablePath}"
        File mutect2Rate = "~{mutect2RatePath}"
    }

    runtime {
        mem: memoryGb + "G"
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/somatic_dna_tools@sha256:f281e73ddf515f3a5db6e766ef12fb2331dafa2b27c88e426764dcd8b4a7e19b"
    }
}

task GetMutect2Rate {
    input {
        Int diskSize = 1
        Int memoryGb = 4
        Array[File] mutect2Logs
        String mutect2RateTablePath
        String mutect2RatePath = "mutect2.rate.txt"
    }

    command {
        python \
        /get_mutect_rate.py \
        ~{sep=" " mutect2Logs} \
        ~{mutect2RateTablePath} \
        > ~{mutect2RatePath}
    }

    output {
        File mutect2RateTable = "~{mutect2RateTablePath}"
        String mutect2Rate = read_float(mutect2RatePath)
    }

    runtime {
        mem: memoryGb + "G"
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/somatic_dna_tools@sha256:f281e73ddf515f3a5db6e766ef12fb2331dafa2b27c88e426764dcd8b4a7e19b"
    }
}

task CompareStrings {
    input {
        String inputString
        File fileInputString
        String fileId
        String resultFilePath = "string_compare_~{fileId}.txt"
        Int diskSize = 1
        Int memoryGb = 1
    }

    command {
        python /compare_strings.py \
        ~{fileInputString} \
        ~{inputString} \
        > ~{resultFilePath}
    }

    output {
        File resultFile = "~{resultFilePath}"
    }

    runtime {
        mem: memoryGb + "G"
        memory : memoryGb + "GB"
        disks: "local-disk " + diskSize + " HDD"
        docker: "gcr.io/nygc-public/somatic_dna_tools@sha256:f281e73ddf515f3a5db6e766ef12fb2331dafa2b27c88e426764dcd8b4a7e19b"
        runtime_minutes: "60"
    }

}

task CreateBlankFile {
    input {
        String fileId
        String blankFilePath = "blank_file_~{fileId}.txt"
        Int diskSize = 1
        Int memoryGb = 1
    }

    command {
        echo "" > ~{blankFilePath}
    }

    output {
        File blankFile = "~{blankFilePath}"
    }

    runtime {
        mem: memoryGb + "G"
        memory : memoryGb + "GB"
        disks: "local-disk " + diskSize + " HDD"
        docker: "gcr.io/nygc-public/genome-utils@sha256:59603ab0aeda571c38811d6d1820d1b546a69fc342120bef75597bfd7905ea1f"
        runtime_minutes: "60"
    }

}

# for wdl version 1.0

task GetIndex {
    input {
        String sampleId
        Array[String] sampleIds
    }

    command {
        python /get_index.py \
        --sample-id ~{sampleId} \
        --sample-ids ~{sep=' ' sampleIds}
    }

    output {
        Int index = read_int(stdout())
    }

    runtime {
        docker: "gcr.io/nygc-public/workflow_utils@sha256:40fa18ac3f9d9f3b9f037ec091cb0c2c26ad6c7cb5c32fb16c1c0cf2a5c9caea"
        runtime_minutes: "60"
    }
}

task SplitVcf {
    input {
        File vcf
        String vcfPath = "~{vcf}"
        String prefix
        Array[String] splitVcfPaths
        String dollarSign = "$"
        Int maxRows = 1000
        Int minSplits = 2
        Int maxSplits = 10
        Int diskSize
        Int memoryGb = 10
    }

    command <<<
        set -e -o pipefail
        rm -f ~{prefix}*.vcf
        
        filename=$( basename  ~{vcfPath} )
        extension="~{dollarSign}{filename##*.}"
        filename="~{dollarSign}{filename%.*}"

        if [[ ~{dollarSign}extension == gz ]]; then
            input_path=~{dollarSign}filename

            gunzip -c \
            ~{vcf} \
            > ~{dollarSign}input_path
        else
            input_path=~{vcf}
        fi
        
        
        python /split_vcf.py \
        --vcf ~{dollarSign}input_path \
        --output-prefix ~{prefix} \
        --max-rows ~{maxRows} \
        --max-splits ~{maxSplits} \
        --min-splits ~{minSplits}
    >>>

    output {
        Array[File?] splitVcfs = splitVcfPaths
    }

    runtime {
        mem: memoryGb + "G"
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-comp-s-fd4e/gcs_htslib_suite@sha256:47eba58683641905a58f31c05605c953dc1d888466e8d434a6ceff47b76df03e"
        runtime_minutes: "1440"
    }
}

task SplitVcfChrom {
    input {
        File vcf
        String prefix
        Array[String] splitVcfPaths
        String chrom
        Int maxRows = 1000
        Int minSplits = 2
        Int maxSplits = 10
        Int diskSize
        Int memoryGb = 10
    }

    command {
        set -e -o pipefail
        rm -f ~{prefix}*.vcf
        
        bgzip -c \
        ~{vcf} \
        > FULL.~{prefix}.vcf.gz
        
        tabix \
        -p vcf \
        FULL.~{prefix}.vcf.gz
        
        bcftools view \
        FULL.~{prefix}.vcf.gz \
        --regions ~{chrom} \
        > CHROM.~{prefix}.vcf

        python /split_vcf.py \
        --vcf CHROM.~{prefix}.vcf \
        --output-prefix ~{prefix} \
        --max-rows ~{maxRows} \
        --max-splits ~{maxSplits} \
        --min-splits ~{minSplits}
    }

    output {
        Array[File?] splitVcfs = splitVcfPaths
    }

    runtime {
        mem: memoryGb + "G"
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-comp-s-fd4e/gcs_htslib_suite@sha256:47eba58683641905a58f31c05605c953dc1d888466e8d434a6ceff47b76df03e"
        runtime_minutes: "1440"
    }
}

task DownloadFile {
    input {
        File uri
        String filePath
        Int diskSize
        Int memoryGb = 10
    }

    command {
        gsutil \
        cp -m \
        ~{uri} \
        .
    }

    output {
        File file = "~{filePath}"
    }

    runtime {
        mem: memoryGb + "G"
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-comp-s-fd4e/gcs_htslib_suite@sha256:47eba58683641905a58f31c05605c953dc1d888466e8d434a6ceff47b76df03e"
        runtime_minutes: "1440"
    }
}

task FilterFile {
    input {
        IndexedVcf referenceVcf
        String filteredReferenceVcfPath
        IndexedReference referenceFa
        Int diskSize
        Int memoryGb = 10
    }

    command {
        set -e -o pipefail
                
        bcftools view \
        ~{referenceVcf.vcf} \
        -f "PASS" \
        |  bcftools annotate \
        -x "^INFO/AF_fin_XX,INFO/AF_nfe_XX,INFO/AF_oth_XX,INFO/AF_sas_XX,INFO/AF_ami_XX,INFO/AF_asj_XX,INFO/AF_mid_XX,INFO/AF_eas_XX,INFO/AF_amr_XX,INFO/AF_afr_XX,INFO/AF_XX,INFO/AF_fin_XY,INFO/AF_nfe_XY,INFO/AF_oth_XY,INFO/AF_sas_XY,INFO/AF_ami_XY,INFO/AF_asj_XY,INFO/AF_mid_XY,INFO/AF_eas_XY,INFO/AF_amr_XY,INFO/AF_afr_XY,INFO/AF_XY,INFO/AF_fin,INFO/AF_nfe,INFO/AF_oth,INFO/AF_sas,INFO/AF_ami,INFO/AF_asj,INFO/AF_mid,INFO/AF_eas,INFO/AF_amr,INFO/AF_afr,INFO/AF,INFO/AF_non_cancer_fin_XX,INFO/AF_non_cancer_nfe_XX,INFO/AF_non_cancer_oth_XX,INFO/AF_non_cancer_sas_XX,INFO/AF_non_cancer_ami_XX,INFO/AF_non_cancer_asj_XX,INFO/AF_non_cancer_mid_XX,INFO/AF_non_cancer_eas_XX,INFO/AF_non_cancer_amr_XX,INFO/AF_non_cancer_afr_XX,INFO/AF_non_cancer_XX,INFO/AF_non_cancer_fin_XY,INFO/AF_non_cancer_nfe_XY,INFO/AF_non_cancer_oth_XY,INFO/AF_non_cancer_sas_XY,INFO/AF_non_cancer_ami_XY,INFO/AF_non_cancer_asj_XY,INFO/AF_non_cancer_mid_XY,INFO/AF_non_cancer_eas_XY,INFO/AF_non_cancer_amr_XY,INFO/AF_non_cancer_afr_XY,INFO/AF_non_cancer_XY,INFO/AF_non_cancer_fin,INFO/AF_non_cancer_nfe,INFO/AF_non_cancer_oth,INFO/AF_non_cancer_sas,INFO/AF_non_cancer_ami,INFO/AF_non_cancer_asj,INFO/AF_non_cancer_mid,INFO/AF_non_cancer_eas,INFO/AF_non_cancer_amr,INFO/AF_non_cancer_afr,INFO/AF_non_cancer,INFO/AN,INFO/AN_non_cancer,INFO/nhomalt,INFO/nhomalt_non_cancer" \
        - \
        |
        bcftools \
        norm \
        -m \
        -any \
        --no-version \
        -f ~{referenceFa.fasta} \
        - \
        | bgzip -c \
        > ~{filteredReferenceVcfPath}

        tabix \
        -p vcf \
        ~{filteredReferenceVcfPath}
    }

    output {
        IndexedVcf filteredReferenceVcf = object {
            vcf : "~{filteredReferenceVcfPath}",
            index : "~{filteredReferenceVcfPath}.tbi"
        }
    }

    runtime {
        mem: memoryGb + "G"
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-comp-s-fd4e/gcs_htslib_suite@sha256:47eba58683641905a58f31c05605c953dc1d888466e8d434a6ceff47b76df03e"
        runtime_minutes: "1440"
    }
}

task FilterExomeFile {
    input {
        IndexedVcf referenceVcf
        String filteredReferenceVcfPath
        IndexedReference referenceFa
        Int diskSize
        Int memoryGb = 10
    }

    command {
        set -e -o pipefail
        
        bcftools view \
        ~{referenceVcf.vcf} \
        -f "PASS" \
        |  bcftools annotate \
        -x "^INFO/AF_fin_female,INFO/AF_nfe_female,INFO/AF_oth_female,INFO/AF_sas_female,INFO/AF_asj_female,INFO/AF_eas_female,INFO/AF_amr_female,INFO/AF_afr_female,INFO/AF_female,INFO/AF_fin_male,INFO/AF_nfe_male,INFO/AF_oth_male,INFO/AF_sas_male,INFO/AF_asj_male,INFO/AF_eas_male,INFO/AF_amr_male,INFO/AF_afr_male,INFO/AF_male,INFO/AF_fin,INFO/AF_nfe,INFO/AF_oth,INFO/AF_sas,INFO/AF_asj,INFO/AF_eas,INFO/AF_amr,INFO/AF_afr,INFO/AF,INFO/non_cancer_AF_fin_female,INFO/non_cancer_AF_nfe_female,INFO/non_cancer_AF_oth_female,INFO/non_cancer_AF_sas_female,INFO/non_cancer_AF_asj_female,INFO/non_cancer_AF_eas_female,INFO/non_cancer_AF_amr_female,INFO/non_cancer_AF_afr_female,INFO/non_cancer_AF_female,INFO/non_cancer_AF_fin_male,INFO/non_cancer_AF_nfe_male,INFO/non_cancer_AF_oth_male,INFO/non_cancer_AF_sas_male,INFO/non_cancer_AF_asj_male,INFO/non_cancer_AF_eas_male,INFO/non_cancer_AF_amr_male,INFO/non_cancer_AF_afr_male,INFO/non_cancer_AF_male,INFO/non_cancer_AF_fin,INFO/non_cancer_AF_nfe,INFO/non_cancer_AF_oth,INFO/non_cancer_AF_sas,INFO/non_cancer_AF_asj,INFO/non_cancer_AF_eas,INFO/non_cancer_AF_amr,INFO/non_cancer_AF_afr,INFO/non_cancer_AF,INFO/AN,INFO/non_cancer_AN,INFO/nhomalt,INFO/non_cancer_nhomalt" \
        - \
        |
        bcftools \
        norm \
        -m \
        -any \
        --no-version \
        -f ~{referenceFa.fasta} \
        - \
        | bgzip -c \
        > ~{filteredReferenceVcfPath}

        tabix \
        -p vcf \
        ~{filteredReferenceVcfPath}
    }

    output {
        IndexedVcf filteredReferenceVcf = object {
            vcf : "~{filteredReferenceVcfPath}",
            index : "~{filteredReferenceVcfPath}.tbi"
        }
    }

    runtime {
        mem: memoryGb + "G"
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-comp-s-fd4e/gcs_htslib_suite@sha256:47eba58683641905a58f31c05605c953dc1d888466e8d434a6ceff47b76df03e"
        runtime_minutes: "1440"
    }
}

task FilterForBiallelicSnps {
    input {
        String gnomadBiallelicPath
        IndexedVcf filteredReferenceVcf

        Int memoryGb = 24
        Int diskSize = (ceil( size(filteredReferenceVcf.vcf, "GB") )  * 2 ) + 10
    }

    command {
        set -e -o pipefail
        
        bcftools view \
        ~{filteredReferenceVcf.vcf} \
        --min-alleles 2 \
        --max-alleles 2 \
        --types snps \
        | bgzip -c \
        > ~{gnomadBiallelicPath}

        tabix \
        -p vcf \
        ~{gnomadBiallelicPath}
    }

    output {
        IndexedVcf gnomadBiallelic = object {
            vcf : "~{gnomadBiallelicPath}",
            index : "~{gnomadBiallelicPath}.tbi"
        }
    }

    runtime {
        mem: memoryGb + "G"
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-comp-s-fd4e/gcs_htslib_suite@sha256:47eba58683641905a58f31c05605c953dc1d888466e8d434a6ceff47b76df03e"
        runtime_minutes: "1440"
    }
}
