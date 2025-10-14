version 1.0

import "../wdl_structs.wdl"

# General tasks

task SliceBam {
    input {
        Int threads = 8
        Int memoryGb = 2
        Int diskSize
        File chromBed
        String sampleId
        String chrom
        String lancetBamSlicePath = "~{sampleId}_~{chrom}.lancet.bam"
        String lancetBaiSlicePath = "~{sampleId}_~{chrom}.lancet.bai"
        Bam finalBam

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

        samtools view \
        --threads ~{threads} \
        -h \
        -b \
        -L ~{chromBed} \
        ~{finalBam.bam} \
        > ~{lancetBamSlicePath}
        
        samtools \
        index \
        -@ ~{threads} \
        ~{lancetBamSlicePath} \
        ~{lancetBaiSlicePath}
    }

    output {
        Bam lancetBamSlice = object {
            bam : lancetBamSlicePath,
            bamIndex : lancetBaiSlicePath
        }
    }

    runtime {
        mem: memoryGb + "G"
        cpus: threads
        cpu : threads
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-comp-s-fd4e/gcs_htslib_suite@sha256:47eba58683641905a58f31c05605c953dc1d888466e8d434a6ceff47b76df03e"
        runtime_minutes: "500"
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

task GetInsertSize {
    input {
        File insertSizeMetrics
    }

    command {
        grep -A 1 \
        "MEDIAN_INSERT_SIZE" \
        ~{insertSizeMetrics} \
        | tail -n 1 \
        | cut -f 1
    }

    output {
        Int insertSize = read_int(stdout())
    }

    runtime {
        docker : "gcr.io/nygc-public/python@sha256:9b7d62026be68c2e91c17fb4e0499454e41ebf498ef345f9ad6e100a67e4b697"
        runtime_minutes: "60"
    }
}

task Gatk4MergeSortVcf {
    input {
        Int diskSize = 10
        Int memoryGb = 8
        String sortedVcfPath
        Array[File] tempChromVcfs
        IndexedReference referenceFa
    }

    Int jvmHeap = memoryGb * 750  # Heap size in Megabytes. mem is in GB. (75% of mem)
    command {
        gatk \
        SortVcf \
        --java-options "-Xmx~{jvmHeap}m -XX:ParallelGCThreads=4" \
        -SD ~{referenceFa.dict} \
        -I ~{sep=" -I " tempChromVcfs} \
        -O ~{sortedVcfPath}
    }

    output {
        IndexedVcf sortedVcf = object {
                vcf : "~{sortedVcfPath}",
                index : "~{sortedVcfPath}.idx"
            }
    }

    parameter_meta {
        sortedVcfPath: {
            description: "Output VCF filename for uncompressed output file. Must not end in .gz. Corresponding index file will end in .idx",
            category: "other"
        }

        tempChromVcfs: {
            description: "Input VCF filenames",
            category: "required"
        }

        referenceFa: {
            description: "Fasta object (the dictionary will be used to order the final VCF)",
            category: "required"
        }
    }

    runtime {
        mem: memoryGb + "G"
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/broadinstitute/gatk4@sha256:b3bde7bc74ab00ddce342bd511a9797007aaf3d22b9cfd7b52f416c893c3774c"
        runtime_minutes: "500"
    }
}

task AddCommandReorderColumnsVcf {
    input {
        Int diskSize
        Int memoryGb
        String normal
        String tumor
        File inVcf
        String orderedVcfPath
        String outVcfPath = sub(sub(basename(inVcf), ".gz$", ""), ".vcf$", "_w_command.vcf")
        File jsonLog
    }

    command {

        python \
        /add_command.py \
        ~{inVcf} \
        ~{outVcfPath} \
        ~{jsonLog}

        python \
        /reorder_vcf.py \
        ~{outVcfPath} \
        ~{orderedVcfPath} \
        ~{normal} ~{tumor}
    }

    output {
        File orderedVcf = "~{orderedVcfPath}"
    }

    runtime {
        mem: memoryGb + "G"
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/somatic_dna_tools@sha256:f281e73ddf515f3a5db6e766ef12fb2331dafa2b27c88e426764dcd8b4a7e19b"
        runtime_minutes: "60"
    }
}

task ReorderVcfColumns {
    input {
        Int diskSize
        Int memoryGb
        String normal
        String tumor
        File rawVcf
        String orderedVcfPath
    }

    command {
        set -e -o pipefail

        if [[ ~{rawVcf} == *.vcf.gz ]]; then
            gunzip \
            -c \
            ~{rawVcf} \
            > ~{rawVcf}.vcf
            runVcf=~{rawVcf}.vcf
        elif [[ ~{rawVcf} == *.vcf.bgz ]]; then
            gunzip \
            -c \
            ~{rawVcf} \
            > ~{rawVcf}.vcf
            runVcf=~{rawVcf}.vcf
        else
            runVcf=~{rawVcf}
        fi

        python \
        /reorder_vcf.py \
        $runVcf \
        ~{orderedVcfPath} \
        ~{normal} ~{tumor}
    }

    output {
        File orderedVcf = "~{orderedVcfPath}"
    }

    runtime {
        mem: memoryGb + "G"
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/somatic_dna_tools@sha256:f281e73ddf515f3a5db6e766ef12fb2331dafa2b27c88e426764dcd8b4a7e19b"
        runtime_minutes: "90"
    }
}

task AddVcfCommand {
    input {
        Int diskSize
        Int memoryGb
        File inVcf
        String outVcfPath = sub(sub(basename(inVcf), ".gz$", ""), ".vcf$", "_w_command.vcf")
        File jsonLog
    }

    command {
        python \
        /add_command.py \
        ~{inVcf} \
        ~{outVcfPath} \
        ~{jsonLog}
    }

    output {
        File outVcf = "~{outVcfPath}"
    }

    runtime {
        mem: memoryGb + "G"
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/somatic_dna_tools@sha256:f281e73ddf515f3a5db6e766ef12fb2331dafa2b27c88e426764dcd8b4a7e19b"
        runtime_minutes: "60"
    }
}

# caller specific tasks

task MantaWgs {
    input {
        Int threads = 8
        Int memoryGb = 4
        Int diskSize
        String pairName
        String intHVmem = "unlimited"
        IndexedReference referenceFa
        Bam normalFinalBam
        IndexedTable callRegions
        Bam tumorFinalBam
    }

    command {
        set -e -o pipefail

        mkdir ~{pairName}.MantaRaw

        configManta.py \
        --normalBam ~{normalFinalBam.bam} \
        --tumorBam ~{tumorFinalBam.bam} \
        --referenceFasta ~{referenceFa.fasta} \
        --callRegions ~{callRegions.table} \
        --runDir ~{pairName}.MantaRaw

        "~{pairName}.MantaRaw/runWorkflow.py" \
        --mode local \
        --job ~{threads} \
        --memGb ~{intHVmem}
    }

    output {
        IndexedVcf candidateSmallIndels = object {
                vcf : "~{pairName}.MantaRaw/results/variants/candidateSmallIndels.vcf.gz",
                index : "~{pairName}.MantaRaw/results/variants/candidateSmallIndels.vcf.gz.tbi"
            }
        IndexedVcf diploidSV = object {
                vcf : "~{pairName}.MantaRaw/results/variants/diploidSV.vcf.gz",
                index : "~{pairName}.MantaRaw/results/variants/diploidSV.vcf.gz.tbi"
            }
        IndexedVcf somaticSV = object {
                vcf : "~{pairName}.MantaRaw/results/variants/somaticSV.vcf.gz",
                index : "~{pairName}.MantaRaw/results/variants/somaticSV.vcf.gz.tbi"
            }
        IndexedVcf candidateSV = object {
                vcf : "~{pairName}.MantaRaw/results/variants/candidateSV.vcf.gz",
                index : "~{pairName}.MantaRaw/results/variants/candidateSV.vcf.gz.tbi"
            }
    }

    runtime {
        mem: memoryGb + "G"
        cpus: threads
        cpu : threads
        disks: "local-disk " + diskSize + " LOCAL"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/manta@sha256:e171112cccf6758693b7a8aab80000fe6121d6969ce83e14a4e15fbc5f2f3662"
        runtime_minutes: "3000"
    }
}

task MantaWgsPon {
    input {
        Int threads = 8
        Int memoryGb = 4
        Int diskSize
        String sampleId
        String intHVmem = "unlimited"
        IndexedReference referenceFa
        IndexedTable callRegions
        Bam tumorFinalBam
    }

    command {
        set -e -o pipefail

        mkdir ~{sampleId}.MantaRaw

        configManta.py \
        --tumorBam ~{tumorFinalBam.bam} \
        --referenceFasta ~{referenceFa.fasta} \
        --callRegions ~{callRegions.table} \
        --runDir ~{sampleId}.MantaRaw

        "~{sampleId}.MantaRaw/runWorkflow.py" \
        --mode local \
        --job ~{threads} \
        --memGb ~{intHVmem}
    }

    output {
        IndexedVcf tumorSV = object {
                vcf : "~{sampleId}.MantaRaw/results/variants/tumorSV.vcf.gz",
                index : "~{sampleId}.MantaRaw/results/variants/tumorSV.vcf.gz.tbi"
            }

    }

    runtime {
        mem: memoryGb + "G"
        cpus: threads
        cpu : threads
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/manta@sha256:e171112cccf6758693b7a8aab80000fe6121d6969ce83e14a4e15fbc5f2f3662"
        runtime_minutes: "3000"
    }
}

task FilterNonpass {
    input {
        Int threads = 4
        Int memoryGb = 8
        Int diskSize
        String pairName
        String outVcfPath = "~{pairName}.manta.filtered.unorder.vcf"
        IndexedReference referenceFa
        File vcf
    }

    Int jvmHeap = memoryGb * 750  # Heap size in Megabytes. mem is in GB. (75% of mem)
    command {
        gatk \
        SelectVariants \
        --java-options "-Xmx~{jvmHeap}m -XX:ParallelGCThreads=4" \
        -R ~{referenceFa.fasta} \
        -V ~{vcf} \
        -O ~{outVcfPath} \
        --exclude-filtered
    }

    output {
        File outVcf = "~{outVcfPath}"
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

task FilterNonpassPon {
    input {
        Int threads = 4
        Int memoryGb = 8
        Int diskSize
        String pairName
        String outVcfPath = "~{pairName}.manta.filtered.unorder.vcf"
        IndexedReference referenceFa
        IndexedVcf vcf
    }

    Int jvmHeap = memoryGb * 750  # Heap size in Megabytes. mem is in GB. (75% of mem)
    command {
        gatk \
        SelectVariants \
        --java-options "-Xmx~{jvmHeap}m -XX:ParallelGCThreads=4" \
        -R ~{referenceFa.fasta} \
        -V ~{vcf.vcf} \
        -O ~{outVcfPath} \
        --exclude-filtered
    }

    output {
        File outVcf = "~{outVcfPath}"
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

task Strelka2 {
    input {
        Int threads
        Int memoryGb
        Int diskSize
        String pairName
        String intHVmem = "unlimited"
        IndexedReference referenceFa
        Bam normalFinalBam
        IndexedTable callRegions
        IndexedVcf candidateSmallIndels
        Bam tumorFinalBam
        File configureStrelkaSomaticWorkflow
    }

    command {
        set -e -o pipefail

        mkdir ~{pairName}.Strelka2Raw

        configureStrelkaSomaticWorkflow.py \
        --normalBam ~{normalFinalBam.bam} \
        --tumorBam ~{tumorFinalBam.bam} \
        --referenceFasta ~{referenceFa.fasta} \
        --callRegions ~{callRegions.table} \
        --indelCandidates ~{candidateSmallIndels.vcf} \
        --config ~{configureStrelkaSomaticWorkflow} \
        --runDir ~{pairName}.Strelka2Raw \

        "~{pairName}.Strelka2Raw/runWorkflow.py" \
        --mode local \
        --job ~{threads} \
        --memGb ~{intHVmem}
    }

    output {
        IndexedVcf strelka2Snvs = object {
                vcf : "~{pairName}.Strelka2Raw/results/variants/somatic.snvs.vcf.gz",
                index : "~{pairName}.Strelka2Raw/results/variants/somatic.snvs.vcf.gz.tbi"
            }
        IndexedVcf strelka2Indels = object {
                vcf : "~{pairName}.Strelka2Raw/results/variants/somatic.indels.vcf.gz",
                index : "~{pairName}.Strelka2Raw/results/variants/somatic.indels.vcf.gz.tbi"
            }
    }

    runtime {
        mem: memoryGb + "G"
        cpus: threads
        cpu : threads
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/strelka@sha256:eb71db1fbe25d67c025251823ae5d5e9dbf5a6861a98298b47ecd722dfa5bd14"
        runtime_minutes: "1000"
    }
}

task Strelka2Exome {
    input {
        Int threads
        Int memoryGb
        Int diskSize
        String pairName
        String intHVmem = "unlimited"
        IndexedReference referenceFa
        Bam normalFinalBam
        IndexedTable callRegions
        Bam tumorFinalBam
        File configureStrelkaSomaticWorkflow
    }

    command {
        set -e -o pipefail

        mkdir ~{pairName}.Strelka2Raw

        configureStrelkaSomaticWorkflow.py \
        --normalBam ~{normalFinalBam.bam} \
        --tumorBam ~{tumorFinalBam.bam} \
        --referenceFasta ~{referenceFa.fasta} \
        --callRegions ~{callRegions.table} \
        --config ~{configureStrelkaSomaticWorkflow} \
        --runDir ~{pairName}.Strelka2Raw \
        --exome

        "~{pairName}.Strelka2Raw/runWorkflow.py" \
        --mode local \
        --job ~{threads} \
        --memGb ~{intHVmem}
    }

    output {
        IndexedVcf strelka2Snvs = object {
                vcf : "~{pairName}.Strelka2Raw/results/variants/somatic.snvs.vcf.gz",
                index : "~{pairName}.Strelka2Raw/results/variants/somatic.snvs.vcf.gz.tbi"
            }
        IndexedVcf strelka2Indels = object {
                vcf : "~{pairName}.Strelka2Raw/results/variants/somatic.indels.vcf.gz",
                index : "~{pairName}.Strelka2Raw/results/variants/somatic.indels.vcf.gz.tbi"
            }
    }

    runtime {
        mem: memoryGb + "G"
        cpus: threads
        cpu : threads
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/strelka@sha256:eb71db1fbe25d67c025251823ae5d5e9dbf5a6861a98298b47ecd722dfa5bd14"
        runtime_minutes: "1000"
    }
}

task SelectVariants {
    input {
        Int threads = 4
        Int memoryGb = 8
        Int diskSize
        String pairName
        String outVcfPath = "~{pairName}.selected.vcf"
        IndexedReference referenceFa
        File intervalListBed
        IndexedVcf vcf
        Int padding = 200
    }

    Int jvmHeap = memoryGb * 750  # Heap size in Megabytes. mem is in GB. (75% of mem)
    command {
        gatk \
        SelectVariants \
        --java-options "-Xmx~{jvmHeap}m -XX:ParallelGCThreads=4" \
        -R ~{referenceFa.fasta} \
        -V ~{vcf.vcf} \
        -O ~{outVcfPath} \
        --intervals ~{intervalListBed} \
        --interval-padding ~{padding}
    }

    output {
        File outVcf = "~{outVcfPath}"
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


task LancetWGSRegional {
    input {
        Int threads
        Int diskSize
        Int memoryGb
        String pairName
        String chrom
        String lancetChromVcfPath = "~{pairName}_~{chrom}.lancet.vcf"
        IndexedReference referenceFa
        File chromBed
        Bam normalFinalBam
        Bam tumorFinalBam
    }

    command {
        lancet \
        --normal ~{normalFinalBam.bam} \
        --tumor ~{tumorFinalBam.bam} \
        --bed ~{chromBed} \
        --ref ~{referenceFa.fasta} \
        --min-k 11 \
        --low-cov 1 \
        --min-phred-fisher 5 \
        --min-strand-bias 1 \
        --min-alt-count-tumor 3 \
        --min-vaf-tumor 0.04 \
        --num-threads ~{threads} \
        > ~{lancetChromVcfPath}
    }

    output {
        File lancetChromVcf = "~{lancetChromVcfPath}"
    }

    runtime {
        mem: memoryGb + "G"
        cpus: threads
        cpu : threads
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/lancet@sha256:25169d34b41de9564e03f02ebcbfb4655cf536449592b0bd58773195f9376e61"
        runtime_minutes: "3000"
    }
}

task LancetExome {
    input {
        Int threads
        Int diskSize
        Int memoryGb
        String pairName
        String chrom
        String lancetChromVcfPath = "~{pairName}_~{chrom}.lancet.vcf"
        IndexedReference referenceFa
        File chromBed
        Bam normalFinalBam
        Bam tumorFinalBam
    }

    command {
        lancet \
        --normal ~{normalFinalBam.bam} \
        --tumor ~{tumorFinalBam.bam} \
        --bed ~{chromBed} \
        --ref ~{referenceFa.fasta} \
        --min-k 11 \
        --low-cov 1 \
        --min-phred-fisher 5 \
        --min-strand-bias 1 \
        --min-alt-count-tumor 3 \
        --min-vaf-tumor 0.04 \
        --num-threads ~{threads} \
        > ~{lancetChromVcfPath}
    }

    output {
        File lancetChromVcf = "~{lancetChromVcfPath}"
    }

    runtime {
        mem: memoryGb + "G"
        cpus: threads
        cpu : threads
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/lancet@sha256:25169d34b41de9564e03f02ebcbfb4655cf536449592b0bd58773195f9376e61"
        runtime_minutes: "3000"
    }
}

task Mutect2Wgs {
    input {
        Int memoryGb = 4
        Int diskSize
        String chrom
        String tumor
        String normal
        String pairName
        IndexedReference referenceFa
        Bam normalFinalBam
        Bam tumorFinalBam
    }

    Int jvmHeap = memoryGb * 750  # Heap size in Megabytes. mem is in GB. (75% of mem)
    String chromFileName = sub(chrom, ":", "-")
    String mutect2ChromRawVcfPath = "~{pairName}_~{chromFileName}.mutect2.raw.vcf"
    String mutect2ChromRawStatsPath = "~{pairName}_~{chromFileName}.mutect2.raw.vcf.stats"
    String mutect2LogPath = "~{pairName}_~{chromFileName}.mutect2.log"

    command {
        #!/bin/bash
        gatk \
        Mutect2 \
        --java-options "-Xmx~{jvmHeap}m -XX:ParallelGCThreads=4" \
        --reference ~{referenceFa.fasta} \
        -L ~{chrom} \
        -I ~{tumorFinalBam.bam} \
        -I ~{normalFinalBam.bam} \
        -tumor ~{tumor} \
        -normal ~{normal} \
        -O ~{mutect2ChromRawVcfPath} \
        2> >(tee -a ~{mutect2LogPath} >&2)
    }

    output {
        File mutect2ChromRawVcf = "~{mutect2ChromRawVcfPath}"
        File mutect2ChromRawStats = "~{mutect2ChromRawStatsPath}"
        File mutect2Log = "~{mutect2LogPath}"
    }

    runtime {
        mem: memoryGb + "G"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/broadinstitute/gatk4@sha256:b3bde7bc74ab00ddce342bd511a9797007aaf3d22b9cfd7b52f416c893c3774c"
        disks: "local-disk " + diskSize + " HDD"
        runtime_minutes: "3000"
    }

    parameter_meta {
        tumorFinalBam: {
            description: "Input tumor BAM",
            category: "required",
            localization_optional: true
        }

        normalFinalBam: {
            description: "Input normal BAM",
            category: "required",
            localization_optional: true
        }

        referenceFa: {
            description: "Fasta object (the dictionary will be used to order the final VCF)",
            category: "required",
            localization_optional: true
        }
    }
}

task Mutect2WgsFewNodes {
    input {
        Int memoryGb = 4
        Int diskSize
        File callerIntervalsBed
        String tumor
        String normal
        String pairName
        IndexedReference referenceFa
        Bam normalFinalBam
        Bam tumorFinalBam
    }

    Int jvmHeap = memoryGb * 750  # Heap size in Megabytes. mem is in GB. (75% of mem)
    String mutect2ChromRawVcfPath = "~{pairName}_FewNodes.mutect2.raw.vcf"
    String mutect2ChromRawStatsPath = "~{pairName}_FewNodes.mutect2.raw.vcf.stats"
    String mutect2LogPath = "~{pairName}_FewNodes.mutect2.log"

    command {
        #!/bin/bash
        gatk \
        Mutect2 \
        --java-options "-Xmx~{jvmHeap}m -XX:ParallelGCThreads=4" \
        --reference ~{referenceFa.fasta} \
        -L ~{callerIntervalsBed} \
        -I ~{tumorFinalBam.bam} \
        -I ~{normalFinalBam.bam} \
        -tumor ~{tumor} \
        -normal ~{normal} \
        -O ~{mutect2ChromRawVcfPath} \
        2> >(tee -a ~{mutect2LogPath} >&2)
    }

    output {
        File mutect2ChromRawVcf = "~{mutect2ChromRawVcfPath}"
        File mutect2ChromRawStats = "~{mutect2ChromRawStatsPath}"
        File mutect2Log = "~{mutect2LogPath}"
    }

    runtime {
        mem: memoryGb + "G"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/broadinstitute/gatk4@sha256:b3bde7bc74ab00ddce342bd511a9797007aaf3d22b9cfd7b52f416c893c3774c"
        disks: "local-disk " + diskSize + " HDD"
        runtime_minutes: "3000"
    }

    parameter_meta {
        tumorFinalBam: {
            description: "Input tumor BAM",
            category: "required",
            localization_optional: true
        }

        normalFinalBam: {
            description: "Input normal BAM",
            category: "required",
            localization_optional: true
        }

        referenceFa: {
            description: "Fasta object (the dictionary will be used to order the final VCF)",
            category: "required",
            localization_optional: true
        }
    }
}

task Mutect2WgsMultiGenotype {
    input {
        Int memoryGb = 4
        Int diskSize
        String splitId
        Array[String] tumors
        String normal
        String participantId
        IndexedReference referenceFa
        IndexedVcf alleles
        Bam normalFinalBam
        Array[Bam] tumorFinalBams
        Array[File] tumorFinalBamList
        String f1r2GzPath = "~{participantId}_~{splitId}.f1r2.tar.gz"
    }

    Int jvmHeap = memoryGb * 750  # Heap size in Megabytes. mem is in GB. (75% of mem)
    String mutect2ChromRawVcfPath = "~{participantId}_~{splitId}.mutect2.raw.vcf"
    String mutect2ChromRawStatsPath = "~{participantId}_~{splitId}.mutect2.raw.vcf.stats"
    String mutect2LogPath = "~{participantId}_~{splitId}.mutect2.log"

    command {
        gatk \
        Mutect2 \
        --java-options "-Xmx~{jvmHeap}m -XX:ParallelGCThreads=4" \
        -alleles ~{alleles.vcf} \
        --force-call-filtered-alleles \
        --f1r2-tar-gz ~{f1r2GzPath} \
        --reference ~{referenceFa.fasta} \
        -L ~{alleles.vcf} \
        -I ~{sep=" -I " tumorFinalBamList} \
        -I ~{normalFinalBam.bam} \
        -normal ~{normal} \
        -O ~{mutect2ChromRawVcfPath} \
        &> ~{mutect2LogPath}
    }

    output {
        File mutect2ChromRawVcf = "~{mutect2ChromRawVcfPath}"
        File mutect2ChromRawStats = "~{mutect2ChromRawStatsPath}"
        File mutect2Log = "~{mutect2LogPath}"
        File f1r2Gz = "~{f1r2GzPath}"
    }

    runtime {
        mem: memoryGb + "G"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/broadinstitute/gatk4@sha256:b3bde7bc74ab00ddce342bd511a9797007aaf3d22b9cfd7b52f416c893c3774c"
        disks: "local-disk " + diskSize + " HDD"
        runtime_minutes: "3000"
    }
    
    parameter_meta {
        tumorFinalBams: {
            description: "Input tumor BAM",
            category: "required",
            localization_optional: true
        }
        
        normalFinalBam: {
            description: "Input normal BAM",
            category: "required",
            localization_optional: true
        }
        
        referenceFa: {
            description: "Fasta object (the dictionary will be used to order the final VCF)",
            category: "required",
            localization_optional: true
        }
    }
}

task LearnReadOrientationModel {
    input {
        Int memoryGb = 16
        Int diskSize = 5
        String pairName
        Array[File] f1r2Gz 
        String outArtifactsPath = "~{pairName}.mutect2.artifacts-priors.tar.gz" 
    }
    Int jvmHeap = memoryGb * 750  # Heap size in Megabytes. mem is in GB. (75% of mem) 
    command {
        gatk \
        LearnReadOrientationModel \
        --java-options "-Xmx~{jvmHeap}m -XX:ParallelGCThreads=4" \
        -I ~{sep=" -I " f1r2Gz} \
        -O ~{outArtifactsPath}

    }

    output {
       File outArtifacts = "~{outArtifactsPath}"
    }

    runtime {
        mem: memoryGb + "G"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/broadinstitute/gatk4@sha256:b3bde7bc74ab00ddce342bd511a9797007aaf3d22b9cfd7b52f416c893c3774c"
        disks: "local-disk " + diskSize + " HDD"
        runtime_minutes: "500"
    }
}

task Mutect2FilterFfpe {
    input {
        Int memoryGb = 4
        Int diskSize = 5
        String pairName
        String mutect2ChromVcfPath = "~{pairName}.filtered.mutect2.vcf"
        IndexedReference referenceFa
        File mutect2ChromRawVcf
        File mutect2ChromRawStats
        File outArtifacts
    }
    
    
   

    Int jvmHeap = memoryGb * 750  # Heap size in Megabytes. mem is in GB. (75% of mem)
    command {
        gatk \
        FilterMutectCalls \
        --java-options "-Xmx~{jvmHeap}m -XX:ParallelGCThreads=4" \
        --reference ~{referenceFa.fasta} \
        -V ~{mutect2ChromRawVcf} \
        --stats ~{mutect2ChromRawStats} \
        --orientation-bias-artifact-priors ~{outArtifacts} \
        -O ~{mutect2ChromVcfPath}
    }
    

    output {
        File mutect2ChromVcf = "~{mutect2ChromVcfPath}"
        
    }

    runtime {
        mem: memoryGb + "G"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/broadinstitute/gatk4@sha256:b3bde7bc74ab00ddce342bd511a9797007aaf3d22b9cfd7b52f416c893c3774c"
        disks: "local-disk " + diskSize + " HDD"
        runtime_minutes: "500"
    }
    
    parameter_meta {

        referenceFa: {
            description: "Fasta object (the dictionary will be used to order the final VCF)",
            category: "required",
            localization_optional: false
        }
    }
}

task Mutect2Exome {
    input {
        Int memoryGb = 4
        Int diskSize
        String chrom
        String tumor
        String normal
        String pairName
        IndexedReference referenceFa
        Bam normalFinalBam
        Bam tumorFinalBam
        Int padding = 200
        File invertedIntervalListBed
    }

    Int jvmHeap = memoryGb * 750  # Heap size in Megabytes. mem is in GB. (75% of mem)
    String chromFileName = sub(chrom, ":", "-")
    String mutect2ChromRawVcfPath = "~{pairName}_~{chromFileName}.mutect2.raw.vcf"
    String mutect2ChromRawStatsPath = "~{pairName}_~{chromFileName}.mutect2.raw.vcf.stats"
    String mutect2LogPath = "~{pairName}_~{chromFileName}.mutect2.log"

    command {
        gatk \
        Mutect2 \
        --java-options "-Xmx~{jvmHeap}m -XX:ParallelGCThreads=4" \
        --reference ~{referenceFa.fasta} \
        --intervals ~{chrom} \
        --interval-padding ~{padding} \
        --exclude-intervals ~{invertedIntervalListBed} \
        -I ~{tumorFinalBam.bam} \
        -I ~{normalFinalBam.bam} \
        -tumor ~{tumor} \
        -normal ~{normal} \
        -O ~{mutect2ChromRawVcfPath} \
        | tee ~{mutect2LogPath}
    }

    output {
        File mutect2ChromRawVcf = "~{mutect2ChromRawVcfPath}"
        File mutect2ChromRawStats = "~{mutect2ChromRawStatsPath}"
        File mutect2Log = "~{mutect2LogPath}"
    }

    runtime {
        mem: memoryGb + "G"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/broadinstitute/gatk4@sha256:b3bde7bc74ab00ddce342bd511a9797007aaf3d22b9cfd7b52f416c893c3774c"
        disks: "local-disk " + diskSize + " HDD"
        runtime_minutes: "3000"
    }

    parameter_meta {
        tumorFinalBam: {
            description: "Input tumor BAM",
            category: "required",
            localization_optional: true
        }

        normalFinalBam: {
            description: "Input normal BAM",
            category: "required",
            localization_optional: true
        }

        referenceFa: {
            description: "Fasta object (the dictionary will be used to order the final VCF)",
            category: "required",
            localization_optional: true
        }
    }
}

task Mutect2ExomeFewNodes {
    input {
        Int memoryGb = 4
        Int diskSize
        File callerIntervalsBed
        String tumor
        String normal
        String pairName
        IndexedReference referenceFa
        Bam normalFinalBam
        Bam tumorFinalBam
        Int padding = 200
        File invertedIntervalListBed
    }

    Int jvmHeap = memoryGb * 750  # Heap size in Megabytes. mem is in GB. (75% of mem)
    String mutect2ChromRawVcfPath = "~{pairName}_FewNodes.mutect2.raw.vcf"
    String mutect2ChromRawStatsPath = "~{pairName}_FewNodes.mutect2.raw.vcf.stats"
    String mutect2LogPath = "~{pairName}_FewNodes.mutect2.log"

    command {
        gatk \
        Mutect2 \
        --java-options "-Xmx~{jvmHeap}m -XX:ParallelGCThreads=4" \
        --reference ~{referenceFa.fasta} \
        --intervals ~{callerIntervalsBed} \
        --interval-padding ~{padding} \
        --exclude-intervals ~{invertedIntervalListBed} \
        -I ~{tumorFinalBam.bam} \
        -I ~{normalFinalBam.bam} \
        -tumor ~{tumor} \
        -normal ~{normal} \
        -O ~{mutect2ChromRawVcfPath} \
        | tee ~{mutect2LogPath}
    }

    output {
        File mutect2ChromRawVcf = "~{mutect2ChromRawVcfPath}"
        File mutect2ChromRawStats = "~{mutect2ChromRawStatsPath}"
        File mutect2Log = "~{mutect2LogPath}"
    }

    runtime {
        mem: memoryGb + "G"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/broadinstitute/gatk4@sha256:b3bde7bc74ab00ddce342bd511a9797007aaf3d22b9cfd7b52f416c893c3774c"
        disks: "local-disk " + diskSize + " HDD"
        runtime_minutes: "3000"
    }

    parameter_meta {
        tumorFinalBam: {
            description: "Input tumor BAM",
            category: "required",
            localization_optional: true
        }

        normalFinalBam: {
            description: "Input normal BAM",
            category: "required",
            localization_optional: true
        }

        referenceFa: {
            description: "Fasta object (the dictionary will be used to order the final VCF)",
            category: "required",
            localization_optional: true
        }
    }
}


task Mutect2WgsPon {
    input {
        Int memoryGb = 4
        Int diskSize
        String chrom
        String tumor
        IndexedReference referenceFa
        Bam tumorFinalBam
    }

    String chromFileName = sub(chrom, ":", "-")
    String mutect2ChromRawStatsPath = "~{tumor}_~{chromFileName}.mutect2.raw.vcf.stats"
    String mutect2ChromRawVcfPath = "~{tumor}_~{chromFileName}.mutect2.raw.vcf"

    command {
        gatk \
        Mutect2 \
        --java-options "-XX:ParallelGCThreads=4" \
        --reference ~{referenceFa.fasta} \
        -L ~{chrom} \
        -I ~{tumorFinalBam.bam} \
        -tumor ~{tumor} \
        -O ~{mutect2ChromRawVcfPath}
    }

    output {
        File mutect2ChromRawVcf = "~{mutect2ChromRawVcfPath}"
        File mutect2ChromRawStats = "~{mutect2ChromRawStatsPath}"
    }

    runtime {
        mem: memoryGb + "G"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/broadinstitute/gatk4@sha256:b3bde7bc74ab00ddce342bd511a9797007aaf3d22b9cfd7b52f416c893c3774c"
        disks: "local-disk " + diskSize + " HDD"
        runtime_minutes: "3000"
    }

    parameter_meta {
        tumorFinalBam: {
            description: "Input tumor BAM",
            category: "required",
            localization_optional: true
        }

        referenceFa: {
            description: "Fasta object (the dictionary will be used to order the final VCF)",
            category: "required",
            localization_optional: true
        }
    }
}

task Mutect2ExomePon {
    input {
        Int memoryGb = 4
        Int diskSize
        String chrom
        String tumor
        IndexedReference referenceFa
        Bam tumorFinalBam
        File invertedIntervalListBed
        Int padding = 200
    }

    String chromFileName = sub(chrom, ":", "-")
    String mutect2ChromRawStatsPath = "~{tumor}_~{chromFileName}.mutect2.raw.vcf.stats"
    String mutect2ChromRawVcfPath = "~{tumor}_~{chromFileName}.mutect2.raw.vcf"

    command {
        gatk \
        Mutect2 \
        --java-options "-XX:ParallelGCThreads=4" \
        --reference ~{referenceFa.fasta} \
        --intervals ~{chrom} \
        --interval-padding ~{padding} \
        --exclude-intervals ~{invertedIntervalListBed} \
        -I ~{tumorFinalBam.bam} \
        -tumor ~{tumor} \
        -O ~{mutect2ChromRawVcfPath}
    }

    output {
        File mutect2ChromRawVcf = "~{mutect2ChromRawVcfPath}"
        File mutect2ChromRawStats = "~{mutect2ChromRawStatsPath}"
    }

    runtime {
        mem: memoryGb + "G"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/broadinstitute/gatk4@sha256:b3bde7bc74ab00ddce342bd511a9797007aaf3d22b9cfd7b52f416c893c3774c"
        disks: "local-disk " + diskSize + " HDD"
        runtime_minutes: "3000"
    }

    parameter_meta {
        tumorFinalBam: {
            description: "Input tumor BAM",
            category: "required",
            localization_optional: true
        }

        referenceFa: {
            description: "Fasta object (the dictionary will be used to order the final VCF)",
            category: "required",
            localization_optional: true
        }
    }
}

task MergeMutectStats {
    input {
        Int memoryGb = 4
        Int diskSize = 5
        Array[File] mutect2ChromRawStats
        String pairName
        String mergedMutect2ChromRawStatsPath = "~{pairName}.mutect2.stats"
    }

    Int jvmHeap = memoryGb * 750  # Heap size in Megabytes. mem is in GB. (75% of mem)
    command {
        gatk \
        MergeMutectStats \
        --java-options "-Xmx~{jvmHeap}m -XX:ParallelGCThreads=4" \
        -stats ~{sep=" -stats " mutect2ChromRawStats} \
        -O ~{mergedMutect2ChromRawStatsPath}
    }

    output {
        File mergedMutect2ChromRawStats = "~{mergedMutect2ChromRawStatsPath}"
    }

    runtime {
        mem: memoryGb + "G"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/broadinstitute/gatk4@sha256:b3bde7bc74ab00ddce342bd511a9797007aaf3d22b9cfd7b52f416c893c3774c"
        disks: "local-disk " + diskSize + " HDD"
        runtime_minutes: "500"
    }
}

task Mutect2Filter {
    input {
        Int memoryGb = 4
        Int diskSize = 5
        String pairName
        String mutect2ChromVcfPath = "~{pairName}.filtered.mutect2.vcf"
        IndexedReference referenceFa
        File mutect2ChromRawVcf
        File mutect2ChromRawStats
    }

    Int jvmHeap = memoryGb * 750  # Heap size in Megabytes. mem is in GB. (75% of mem)
    command {
        gatk \
        FilterMutectCalls \
        --java-options "-Xmx~{jvmHeap}m -XX:ParallelGCThreads=4" \
        --reference ~{referenceFa.fasta} \
        -V ~{mutect2ChromRawVcf} \
        --stats ~{mutect2ChromRawStats} \
        -O ~{mutect2ChromVcfPath}
    }

    output {
        File mutect2ChromVcf = "~{mutect2ChromVcfPath}"
    }

    runtime {
        mem: memoryGb + "G"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/broadinstitute/gatk4@sha256:b3bde7bc74ab00ddce342bd511a9797007aaf3d22b9cfd7b52f416c893c3774c"
        disks: "local-disk " + diskSize + " HDD"
        runtime_minutes: "500"
    }

    parameter_meta {

        referenceFa: {
            description: "Fasta object (the dictionary will be used to order the final VCF)",
            category: "required",
            localization_optional: false
        }
    }
}

task SvabaWgs {
    input {
        Int threads
        Int memoryGb = 16
        String pairName
        IndexedTable callRegions
        BwaReference bwaReference
        Bam normalFinalBam
        File dbsnpIndels
        Bam tumorFinalBam
        Int diskSize
        Int preemptible = 3
    }

    command {
        svaba \
        run \
        -t ~{tumorFinalBam.bam} \
        -n ~{normalFinalBam.bam} \
        -p ~{threads} \
        --region ~{callRegions.table} \
        -D ~{dbsnpIndels} \
        -a ~{pairName} \
        -G ~{bwaReference.fasta} \
        -z on
    }

    output {
        File svabaRawGermlineIndel = "~{pairName}.svaba.germline.indel.vcf.gz"
        File svabaRawGermlineSv = "~{pairName}.svaba.germline.sv.vcf.gz"
        File svabaIndelGz = "~{pairName}.svaba.somatic.indel.vcf.gz"
        File svabaGz = "~{pairName}.svaba.somatic.sv.vcf.gz"
    }

    runtime {
        mem: memoryGb + "G"
        cpus: threads
        cpu : threads
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/svaba@sha256:48f6bd86e933ca88fd74d8effc66e93eee5b40945ee37612b80d7edaadc567f3"
        preemptible: preemptible
        runtime_minutes: "4000"
    }
}

task PopulateCache {
    input {
        Int memoryGb = 16
        BwaReference bwaReference
        String refCacheDirPath = "ref_cache"
        String refCachePath = "ref_cache.tar.gz"
        Int diskSize = 10
    }

    command {
        set -e -o pipefail

        /samtools-1.4.1/misc/seq_cache_populate.pl \
        -root ~{refCacheDirPath} \
        ~{bwaReference.fasta}

        tar -czvf \
        ~{refCachePath} \
        ~{refCacheDirPath}
    }

    output {
        File refCache = "~{refCachePath}"
    }

    runtime {
        mem: memoryGb + "G"
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/samtools@sha256:e1149e965e8379f4a75b120d832b84e87dbb97bd5510ed581113400f768e5940"
        runtime_minutes: "500"
    }
}

task SvabaIndex {
    input {
        Int memoryGb = 16
        BwaReference bwaReference
        String svabaIndexedReferencePath = basename(bwaReference.fasta)
        Int diskSize = 16
    }

    command {
        set -e -o pipefail

        cp ~{bwaReference.fasta} \
        ~{svabaIndexedReferencePath}

        /svaba/SeqLib/bwa/bwa \
        index \
        ~{svabaIndexedReferencePath}
    }

    output {
        BwaReference svabaIndexedReference = object {
                fasta : "~{svabaIndexedReferencePath}",
                sa : "~{svabaIndexedReferencePath}.sa",
                pac : "~{svabaIndexedReferencePath}.pac",
                bwt : "~{svabaIndexedReferencePath}.bwt",
                ann : "~{svabaIndexedReferencePath}.ann",
                amb : "~{svabaIndexedReferencePath}.amb"
            }
    }

    runtime {
        mem: memoryGb + "G"
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/svaba@sha256:48f6bd86e933ca88fd74d8effc66e93eee5b40945ee37612b80d7edaadc567f3"
        runtime_minutes: "500"
    }
}

task SvabaWgsPon {
    input {
        Int threads = 4
        Int memoryGb = 16
        String sampleId
        IndexedTable callRegions
        BwaReference svabaIndexedReference
        File dbsnpIndels
        File refCache
        Bam tumorFinalBam
        String refCacheDirPath = sub(basename(refCache), ".tar.gz$", "")
        Int diskSize
        Int verbose = 0
        Int preemptible = 3
    }

    command {
        set -e -o pipefail

        tar -xzf \
        ~{refCache}

        export REF_PATH=./~{refCacheDirPath}/%2s/%2s/%s
        export REF_CACHE=./~{refCacheDirPath}/%2s/%2s/%s

        svaba \
        run \
        --verbose ~{verbose} \
        -t ~{tumorFinalBam.bam} \
        -p ~{threads} \
        -L 100000 \
        --region ~{callRegions.table} \
        -D ~{dbsnpIndels} \
        -a ~{sampleId} \
        -G ~{svabaIndexedReference.fasta} \
        -z on
    }

    output {
        File svabaIndelGz = "~{sampleId}.svaba.indel.vcf.gz"
        File svabaGz = "~{sampleId}.svaba.sv.vcf.gz"
    }

    runtime {
        mem: memoryGb + "G"
        cpus: threads
        cpu : threads
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/svaba@sha256:48f6bd86e933ca88fd74d8effc66e93eee5b40945ee37612b80d7edaadc567f3"
        preemptible: preemptible
        runtime_minutes: "500"
    }
}

task SvabaWgsPonNoL {
    input {
        Int threads = 4
        Int memoryGb = 16
        String sampleId
        IndexedTable callRegions
        BwaReference svabaIndexedReference
        File dbsnpIndels
        File refCache
        Bam tumorFinalBam
        String refCacheDirPath = sub(basename(refCache), ".tar.gz$", "")
        Int diskSize
        Int verbose = 0
        Int preemptible = 3
    }

    command {
        set -e -o pipefail

        tar -xzf \
        ~{refCache}

        export REF_PATH=./~{refCacheDirPath}/%2s/%2s/%s
        export REF_CACHE=./~{refCacheDirPath}/%2s/%2s/%s

        svaba \
        run \
        --verbose ~{verbose} \
        -t ~{tumorFinalBam.bam} \
        -p ~{threads} \
        --region ~{callRegions.table} \
        -D ~{dbsnpIndels} \
        -a ~{sampleId} \
        -G ~{svabaIndexedReference.fasta} \
        -z on
    }

    output {
        File svabaIndelGz = "~{sampleId}.svaba.indel.vcf.gz"
        File svabaGz = "~{sampleId}.svaba.sv.vcf.gz"
    }

    runtime {
        mem: memoryGb + "G"
        cpus: threads
        cpu : threads
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/svaba@sha256:48f6bd86e933ca88fd74d8effc66e93eee5b40945ee37612b80d7edaadc567f3"
        preemptible: preemptible
        runtime_minutes: "500"
    }
}

task ReheaderVcf {
    input {
        Array[String] sampleIds
        String outVcfPath
        IndexedReference referenceFa
        File inVcf
        String outType = "v"
        Int memoryGb = 16
        Int diskSize = (ceil( size(inVcf, "GB") )  * 2 ) + 1
    }

    command {
        set -e -o pipefail

        for id in ~{sep=" " sampleIds} ; do
            echo $id >> sample_list.txt
        done

        bcftools \
        reheader \
        --samples sample_list.txt \
        ~{inVcf} \
        | bcftools \
        convert \
        -O ~{outType} \
        --output ~{outVcfPath} \
        -
    }

    output {
        File reheaderedVcf = "~{outVcfPath}"
    }

    runtime {
        mem: memoryGb + "G"
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/bcftools@sha256:d9b84254b8cc29fcae76b728a5a9a9a0ef0662ee52893d8d82446142876fb400"
        runtime_minutes: "60"
    }
}

task UniqReads {
    input {
        Int memoryGb
        String sampleId
        String tempPrefix = "~{sampleId}_"
        Array[String] tempSeqsPaths
        Bam finalBam
        Int diskSize
    }

    command {
        /samtools-0.1.7a_getUnique-0.1.3/samtools \
        view \
        -U "BWA,~{tempPrefix},N,N" \
        ~{finalBam.bam}

    }

    output {
        Array[File] tempSeqs = tempSeqsPaths
    }

    runtime {
        mem: memoryGb + "G"
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/bicseq2@sha256:3d110b672df0385f761fb64fcf63e98685e6d810c5560043efed5ab94961f3a9"
        runtime_minutes: "500"
    }
}


task Bicseq2Norm {
    input {
        Int memoryGb
        Int readLength
        Int medianInsertSize
        String GCvsRDPath = "~{sampleId}.GCvsRD.pdf"
        String paramsPath = "~{sampleId}.params.out"
        String sampleId

        Array[File] tempSeqs
        Array[String] tempNormPaths

        File bicseq2ConfigFile
        String configFilePath = "~{sampleId}.bicseq2.config"

        Array[File] uniqCoordsFiles
        Array[File] chromFastasFiles
        Int diskSize = 70
    }

    command {
        set -e -o pipefail

        python3 \
        /bicseq2_config_writer.py \
        --fa-files ~{sep=" " chromFastasFiles} \
        --mappability-files ~{sep=" " uniqCoordsFiles} \
        --temp-seqs ~{sep=" " tempSeqs} \
        --temp-norm-paths ~{sep=" " tempNormPaths} \
        --norm-bicseq2-config ~{bicseq2ConfigFile} \
        --sample-id ~{sampleId} \
        --out-file ~{configFilePath}

        mkdir -p ~{sampleId}

        /NBICseq-norm_v0.2.4/NBICseq-norm.pl \
        -l=~{readLength} \
        -s=~{medianInsertSize} \
        -fig=~{GCvsRDPath} \
        -tmp=~{sampleId} \
        ~{configFilePath} \
        ~{paramsPath}
    }

    output {
        Array[File] tempNorm = tempNormPaths
        File params = paramsPath
    }

    runtime {
        mem: memoryGb + "G"
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/bicseq2@sha256:3d110b672df0385f761fb64fcf63e98685e6d810c5560043efed5ab94961f3a9"
        runtime_minutes: "500"
    }
}

task Bicseq2NormReader {
    input {
        Int memoryGb
        Int readLength
        File insertSizeMetrics
        String GCvsRDPath = "~{sampleId}.GCvsRD.pdf"
        String paramsPath = "~{sampleId}.params.out"
        String sampleId

        Array[File] tempSeqs
        Array[String] tempNormPaths

        File bicseq2ConfigFile
        String configFilePath = "~{sampleId}.bicseq2.config"

        Array[File] uniqCoordsFiles
        Array[File] chromFastasFiles
        String dollarSign = "$"
        Int diskSize = 70
    }

    command {
        set -e -o pipefail
        
        medianInsertSize=~{dollarSign}(grep -A 1 "MEDIAN_INSERT_SIZE" ~{insertSizeMetrics} | tail -n 1 | cut -f 1 )

        python3 \
        /bicseq2_config_writer.py \
        --fa-files ~{sep=" " chromFastasFiles} \
        --mappability-files ~{sep=" " uniqCoordsFiles} \
        --temp-seqs ~{sep=" " tempSeqs} \
        --temp-norm-paths ~{sep=" " tempNormPaths} \
        --norm-bicseq2-config ~{bicseq2ConfigFile} \
        --sample-id ~{sampleId} \
        --out-file ~{configFilePath}

        mkdir -p ~{sampleId}

        /NBICseq-norm_v0.2.4/NBICseq-norm.pl \
        -l=~{readLength} \
        -s=~{dollarSign}medianInsertSize \
        -fig=~{GCvsRDPath} \
        -tmp=~{sampleId} \
        ~{configFilePath} \
        ~{paramsPath}
    }

    output {
        Array[File] tempNorm = tempNormPaths
        File params = paramsPath
    }

    runtime {
        mem: memoryGb + "G"
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/bicseq2@sha256:3d110b672df0385f761fb64fcf63e98685e6d810c5560043efed5ab94961f3a9"
        runtime_minutes: "500"
    }
}


task Bicseq2Wgs {
    input {
        Int memoryGb
        String pairName
        String bicseq2PngPath = "~{pairName}.bicseq2.png"
        String bicseq2Path = "~{pairName}.bicseq2.txt"
        Array[File] tempTumorNorms
        Array[File] tempNormalNorms
        File bicseq2SegConfigFile
        String segConfigFilePath = "~{pairName}.bicseq2.seg.config"
        Int lambda = 4
        Int diskSize = 10
    }

    command {
        set -e -o pipefail

        mkdir -p ~{pairName}

        python3 \
        /bicseq2_seg_config_writer.py \
        --tumor-norms ~{sep=" " tempTumorNorms} \
        --normal-norms ~{sep=" " tempNormalNorms} \
        --seg-bicseq2-config ~{bicseq2SegConfigFile} \
        --out-file ~{segConfigFilePath} \
        --pair-id ~{pairName}

        perl /NBICseq-seg_v0.7.2/NBICseq-seg.pl \
        --control \
        --tmp ~{pairName} \
        --fig=~{bicseq2PngPath} \
        --title=~{pairName} \
        --lambda=4 \
        ~{segConfigFilePath} \
        ~{bicseq2Path}
    }

    output {
        File bicseq2Png = "~{pairName}.bicseq2.png"
        File bicseq2 = "~{pairName}.bicseq2.txt"
    }

    runtime {
        mem: memoryGb + "G"
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/bicseq2@sha256:3d110b672df0385f761fb64fcf63e98685e6d810c5560043efed5ab94961f3a9"
        runtime_minutes: "60"
    }
}


task GridssPreprocess {
    input {
        Int threads
        Int memoryGb = 16
        Int diskSize = 700
        Bam finalBam

        BwaReference bwaReference
        Array[File] gridssAdditionalReference
        Int timeLimit = 1000
    }

    String bamBase = basename(finalBam.bam)
    String svBamPath = bamBase + ".gridss.working/" + bamBase + ".sv.bam"
    String cigarMetricsPath = bamBase + ".gridss.working/" + bamBase + ".cigar_metrics"
    String idsvMetricsPath = bamBase + ".gridss.working/" + bamBase + ".idsv_metrics"
    String tagMetricsPath = bamBase + ".gridss.working/" + bamBase + ".tag_metrics"
    String mapqMetricsPath = bamBase + ".gridss.working/" + bamBase + ".mapq_metrics"
    String insertSizeMetricsPath = bamBase + ".gridss.working/" + bamBase + ".insert_size_metrics"

    command {
        set -e -o pipefail

        # mv gridss fasta refs into position in the fasta dir
        fasta_dir=$( dirname ~{bwaReference.fasta} )
        mv ~{sep=" " gridssAdditionalReference} $fasta_dir

        working=$( pwd )

        timeout ~{timeLimit}m \
        gridss \
        --steps preprocess \
        --reference ~{bwaReference.fasta} \
        --jar /opt/gridss/gridss-2.13.2-gridss-jar-with-dependencies.jar \
        --threads ~{threads} \
        --workingdir $working \
        --picardoptions VALIDATION_STRINGENCY=LENIENT \
        ~{finalBam.bam} \
        && cat *.log
    }

    output {
        Bam svBam = object {
            bam : svBamPath,
            bamIndex : sub(svBamPath, ".bam$", ".bam.csi")
        }
        File cigarMetrics = "~{cigarMetricsPath}"
        File idsvMetrics = "~{idsvMetricsPath}"
        File tagMetrics = "~{tagMetricsPath}"
        File mapqMetrics = "~{mapqMetricsPath}"
        File insertSizeMetrics = "~{insertSizeMetricsPath}"
    }

    runtime {
        mem: memoryGb + "G"
        cpus: threads
        disks: "local-disk " + diskSize + " HDD"
        cpu : threads
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/gridss@sha256:bf4a0367accd33911a5658491e42e5ee52440cd196ea96e538511e88a71dee1f"
        runtime_minutes: "1000"
    }
}

task GridssAssembleChunk {
    input {
        Int threads
        Int memoryGb = 48
        Int diskSize = 700
        String pairName

        String gridssassemblyBamPath = "~{pairName}.gridssassembly.bam"
        Int jobIndex
        Int assembleChunks
        BwaReference bwaReference
        Array[File] gridssAdditionalReference
        Bam tumorFinalBam
        Bam normalFinalBam

        Bam normalSvBam
        File normalCigarMetrics
        File normalIdsvMetrics
        File normalTagMetrics
        File normalMapqMetrics
        File normalInsertSizeMetrics

        Bam tumorSvBam
        File tumorCigarMetrics
        File tumorIdsvMetrics
        File tumorTagMetrics
        File tumorMapqMetrics
        File tumorInsertSizeMetrics
    }

    command {
        set -e -o pipefail

        # mv gridss fasta refs into position in the fasta dir
        fasta_dir=$( dirname ~{bwaReference.fasta} )
        mv ~{sep=" " gridssAdditionalReference} $fasta_dir


        # reposition unnamed input
        bash gridss_arrange.sh \
        --tumorFinalBam ~{tumorFinalBam.bam} \
        --normalFinalBam ~{normalFinalBam.bam} \
        --tumorSvBam ~{tumorSvBam.bam} \
        --normalSvBam ~{normalSvBam.bam}

        working=$( pwd )

        bash gridss \
        --steps assemble \
        --reference ~{bwaReference.fasta} \
        --jar /opt/gridss/gridss-2.13.2-gridss-jar-with-dependencies.jar \
        --threads ~{threads} \
        --workingdir $working \
        --assembly ~{gridssassemblyBamPath} \
        --picardoptions VALIDATION_STRINGENCY=LENIENT \
        --jobindex ~{jobIndex} \
        --jobnodes ~{assembleChunks} \
        ~{normalFinalBam.bam} ~{tumorFinalBam.bam} \
        && cat *.log
    }

    output {
        File downsampled = "~{pairName}.gridssassembly.bam.gridss.working/~{pairName}.gridssassembly.bam.downsampled_~{jobIndex}.bed"
        File excluded = "~{pairName}.gridssassembly.bam.gridss.working/~{pairName}.gridssassembly.bam.excluded_~{jobIndex}.bed"
        File subsetCalled = "~{pairName}.gridssassembly.bam.gridss.working/~{pairName}.gridssassembly.bam.subsetCalled_~{jobIndex}.bed"
    }

    runtime {
        mem: memoryGb + "G"
        cpus: threads
        disks: "local-disk " + diskSize + " LOCAL"
        cpu : threads
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/gridss@sha256:bf4a0367accd33911a5658491e42e5ee52440cd196ea96e538511e88a71dee1f"
        runtime_minutes: "1000"
    }
}


task GridssAssemble {
    input {
        Int threads
        Int memoryGb = 48
        Int diskSize = 700
        String pairName
        String gridssassemblyBamPath = "~{pairName}.gridssassembly.bam"
        String gridssassemblySvBamPath = "~{pairName}.gridssassembly.bam.gridss.working/~{pairName}.gridssassembly.bam.sv.bam"
        BwaReference bwaReference
        Array[File] gridssAdditionalReference
        Bam tumorFinalBam
        Bam normalFinalBam
        Array[File] downsampled
        Array[File] excluded
        Array[File] subsetCalled

        # from earlier
        Bam normalSvBam
        File normalCigarMetrics
        File normalIdsvMetrics
        File normalTagMetrics
        File normalMapqMetrics
        File normalInsertSizeMetrics
        Bam tumorSvBam
        File tumorCigarMetrics
        File tumorIdsvMetrics
        File tumorTagMetrics
        File tumorMapqMetrics
        File tumorInsertSizeMetrics
    }

    command {
        set -e -o pipefail

        # mv gridss fasta refs into position in the fasta dir
        fasta_dir=$( dirname ~{bwaReference.fasta} )
        mv ~{sep=" " gridssAdditionalReference} $fasta_dir

        # link preprocess results
        # reposition unnamed input
        bash gridss_arrange.sh \
        --tumorFinalBam ~{tumorFinalBam.bam} \
        --normalFinalBam ~{normalFinalBam.bam} \
        --tumorSvBam ~{tumorSvBam.bam} \
        --normalSvBam ~{normalSvBam.bam}

        # link chunk assembly results
        sub_dir=~{gridssassemblyBamPath}.gridss.working/
        mkdir -p $sub_dir
        copied_sub_dir=$( dirname ~{downsampled[0]} )
        ln -s $copied_sub_dir/* $sub_dir

        working=$( pwd )

        bash gridss \
        --steps assemble \
        --reference ~{bwaReference.fasta} \
        --jar /opt/gridss/gridss-2.13.2-gridss-jar-with-dependencies.jar \
        --threads ~{threads} \
        --workingdir $working \
        --assembly ~{gridssassemblyBamPath} \
        --picardoptions VALIDATION_STRINGENCY=LENIENT \
        ~{normalFinalBam.bam} ~{tumorFinalBam.bam} \
        && cat *.log
    }

    output {
        File gridssassemblyBam = gridssassemblyBamPath
        Bam gridssassemblySvBam = object {
            bam : gridssassemblySvBamPath,
            bamIndex : sub(gridssassemblySvBamPath, ".bam$", ".bam.bai")
        }
    }

    runtime {
        mem: memoryGb + "G"
        cpus: threads
        disks: "local-disk " + diskSize + " LOCAL"
        cpu : threads
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/gridss@sha256:bf4a0367accd33911a5658491e42e5ee52440cd196ea96e538511e88a71dee1f"
        runtime_minutes: "3000"
    }
}

task GridssCalling {
    input {
        Int threads
        Int memoryGb = 48
        Int diskSize = 700
        String pairName
        String gridssassemblyBamPath = "~{pairName}.gridssassembly.bam"
        String gridssUnfilteredVcfPath = "~{pairName}.sv.gridss.unfiltered.vcf"
        BwaReference bwaReference
        Array[File] gridssAdditionalReference
        Array[File] downsampled
        Array[File] excluded
        Array[File] subsetCalled
        File gridssassemblyBam
        Bam gridssassemblySvBam
        Bam tumorFinalBam
        Bam normalFinalBam

        # required from Preprocessing
        Bam normalSvBam
        File normalCigarMetrics
        File normalIdsvMetrics
        File normalTagMetrics
        File normalMapqMetrics
        File normalInsertSizeMetrics

        Bam tumorSvBam
        File tumorCigarMetrics
        File tumorIdsvMetrics
        File tumorTagMetrics
        File tumorMapqMetrics
        File tumorInsertSizeMetrics

        Int timeLimit = 1000
    }

    command {
        set -e -o pipefail

        # mv gridss fasta refs into position in the fasta dir
        fasta_dir=$( dirname ~{bwaReference.fasta} )
        mv ~{sep=" " gridssAdditionalReference} $fasta_dir

        # link assembly results
        sub_dir=~{gridssassemblyBamPath}.gridss.working/
        mkdir -p $sub_dir
        copied_sub_dir=$( dirname ~{gridssassemblySvBam.bam} )
        ln -s $copied_sub_dir/* $sub_dir

        ln -s \
        ~{gridssassemblyBam} \
        ~{gridssassemblyBamPath}


        # link chunk assembly results
        sub_dir=~{gridssassemblyBamPath}.gridss.working/
        mkdir -p $sub_dir
        copied_sub_dir=$( dirname ~{downsampled[0]} )
        ln -s $copied_sub_dir/* $sub_dir

        # link preprocess results
        # reposition unnamed input
        bash gridss_arrange.sh \
        --tumorFinalBam ~{tumorFinalBam.bam} \
        --normalFinalBam ~{normalFinalBam.bam} \
        --tumorSvBam ~{tumorSvBam.bam} \
        --normalSvBam ~{normalSvBam.bam}

        working=$( pwd )

        timeout ~{timeLimit}m \
        bash gridss \
        --steps call \
        --reference ~{bwaReference.fasta} \
        --jar /opt/gridss/gridss-2.13.2-gridss-jar-with-dependencies.jar \
        --threads ~{threads} \
        --workingdir $working \
        --assembly ~{gridssassemblyBamPath} \
        --output ~{gridssUnfilteredVcfPath} \
        --picardoptions VALIDATION_STRINGENCY=LENIENT \
        ~{normalFinalBam.bam} ~{tumorFinalBam.bam} \
        && cat *.log
    }

    output {
        File gridssUnfilteredVcf = "~{pairName}.sv.gridss.unfiltered.vcf"
    }

    runtime {
        mem: memoryGb + "G"
        cpus: threads
        disks: "local-disk " + diskSize + " HDD"
        cpu : threads
        memory : memoryGb + "GB"
        docker: "gcr.io/nygc-public/gridss@sha256:bf4a0367accd33911a5658491e42e5ee52440cd196ea96e538511e88a71dee1f"
        runtime_minutes: "1000"
    }
}

task FilterNonChroms {
    input {
        Int diskSize
        Int memoryGb
        File gridssUnfilteredVcf
        String pairName
        String gridssUnfilteredVcfChromsPath = "~{pairName}.sv.gridss.unfiltered.chroms.vcf"
        Array[String]+ listOfChroms
    }

    command {
        python \
        /vcf_filter.py \
        --vcf-file ~{gridssUnfilteredVcf} \
        --output ~{gridssUnfilteredVcfChromsPath} \
        --chroms ~{sep=" " listOfChroms}
    }

    output {
        File gridssUnfilteredVcfChroms = "~{gridssUnfilteredVcfChromsPath}"
    }

    runtime {
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/somatic_dna_tools@sha256:f281e73ddf515f3a5db6e766ef12fb2331dafa2b27c88e426764dcd8b4a7e19b"
        runtime_minutes: "60"
    }
}

task GridssFilter {
    input {
        Int threads
        Int memoryGb = 16
        Int diskSize = 4
        String bsGenome
        String pairName
        File ponTarGz
        String gridssVcfPath = "~{pairName}.sv.gridss.vcf"
        String gridssVcfPathOut = "~{pairName}.sv.gridss.vcf.bgz"
        String tumourordinal = 2
        File gridssUnfilteredVcf
    }

    command {
        set -e -o pipefail

        tar -zxvf ~{ponTarGz}

        working=$( pwd )

        Rscript /opt/gridss/scripts/gridss_somatic_filter \
        --ref ~{bsGenome} \
        --input ~{gridssUnfilteredVcf} \
        --output ~{gridssVcfPath} \
        --tumourordinal ~{tumourordinal} \
        --plotdir $working \
        --scriptdir /opt/gridss/scripts/ \
        --configdir /opt/gridss/scripts/ \
        --pondir pon/
    }

    output {
        IndexedVcf gridssVcf = object {
                vcf : "~{gridssVcfPathOut}",
                index : "~{gridssVcfPathOut}.tbi"
            }
    }

    runtime {
        mem: memoryGb + "G"
        cpus: threads
        disks: "local-disk " + diskSize + " HDD"
        cpu : threads
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/gridss_filter@sha256:08a1a6e3af74ae8f2a4b84c41f707a05c7271bddcd011d1e63edb0ec88366e11"
        runtime_minutes: "500"
    }
}
