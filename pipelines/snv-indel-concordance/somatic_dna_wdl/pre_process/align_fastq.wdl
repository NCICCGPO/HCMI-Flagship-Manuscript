version 1.0

import "../wdl_structs.wdl"

task Skewer {
    input {
        # command
        Fastqs fastqs
        File adaptersFa
        # replace .fastq.gz or .fq.gz with "-trimmed-pair2.fastq.gz"
        String fastqPrefix = sub(basename(fastqs.fastqR1), ".\\.f.*q\\.gz$", "")
        String fastqOutR1Path = sub(basename(fastqs.fastqR1), ".\\.f.*q\\.gz$", "-trimmed-pair1.fastq.gz")
        String fastqOutR2Path = sub(basename(fastqs.fastqR1), ".\\.f.*q\\.gz$", "-trimmed-pair2.fastq.gz")
        # resources
        Int threads = 16
        Int memoryGb = 24
        Int diskSize
        Int timeLimit = 500
    }

    command {
        timeout ~{timeLimit}m \
        skewer \
        --compress \
        -x ~{adaptersFa} \
        -f "sanger" \
        -t ~{threads} \
        -m "pe" \
        -o ~{fastqPrefix} \
        ~{fastqs.fastqR1} \
        ~{fastqs.fastqR2}
    }

    output {
        Fastqs skewerFastqs = object {
            fastqR1 : fastqOutR1Path,
            fastqR2 : fastqOutR2Path,
            clientSampleId: fastqs.clientSampleId,
            limsLibraryName: fastqs.limsLibraryName,
            readGroupId: fastqs.readGroupId,
            readGroupPlatformUnit: fastqs.readGroupPlatformUnit
        }
    }

    runtime {
        mem: memoryGb + "G"
        cpus: threads
        cpu : threads
        memory : memoryGb + " GB"
        docker : "gcr.io/nygc-public/skewer@sha256:31cb01e3deab46a7a7e9e93e789ab68e66d5b660b7d0b8967c096da6fd38c5e9"
        disks: "local-disk " + diskSize + " HDD"
        runtime_minutes: "500"
    }
}

task AlignBwaMem2 {
    input {
        # command
        Fastqs fastqsAlign
        BwaMem2Reference bwamem2Reference
        String laneBamPath = "~{fastqsAlign.readGroupId}.readgroup.bam"
        # resources
        Int memoryGb
        Int threads = 4
        Int bwamem2Threads = 32
        Int totalThreads = threads + bwamem2Threads
        Int diskSize

        # Values used in RG tags. These are overridden for external fastqs or if we start using
        # other sequencing platforms.
        String platform = "illumina"
        String machineType = "NovaSeq"
        String center = "NYGC"
    }
    String libraryName = select_first([fastqsAlign.limsLibraryName, fastqsAlign.clientSampleId])
    command {
        set -e -o pipefail
        bwa-mem2 mem \
        -Y \
        -K 100000000 \
        -t ~{bwamem2Threads} \
        -R '@RG\tID:~{fastqsAlign.readGroupId}\tPL:~{platform}\tPM:~{machineType}\tLB:~{libraryName}\tDS:hg38\tSM:~{fastqsAlign.clientSampleId}\tCN:~{center}\tPU:${fastqsAlign.readGroupPlatformUnit}' \
        ~{bwamem2Reference.fasta} \
        ~{fastqsAlign.fastqR1} \
        ~{fastqsAlign.fastqR2} \
        | samtools view \
        -@ ~{threads} \
        -Shb \
        -o ~{laneBamPath} \
        -
    }
    output {
        File laneBam = laneBamPath
    }
    runtime {
        mem: memoryGb + "G"
        cpus: totalThreads
        cpu : totalThreads
        memory : memoryGb + " GB"
        disks: "local-disk " + diskSize + " LOCAL"
        docker : "gcr.io/nygc-public/bwamem2@sha256:d46d1433562b8a4e56486298ef522f02f17cb7bd079962b246245136cde15743"
        runtime_minutes: "500"
    }
    meta {
        tool : "bwa-mem"
        version : "2-2.2.1_x64-linux"
    }
}

task AlignMinimap2 {
    input {
        # command
        Fastqs fastqsAlign
        BwaReference bwaReference
        String laneBamPath = "~{fastqsAlign.readGroupId}.readgroup.bam"
        # resources
        Int memoryGb
        Int threads
        Int minimapThreads
        Int totalThreads = 80
        Int diskSize

        # Values used in RG tags. These are overridden for external fastqs or if we start using
        # other sequencing platforms.
        String platform = "illumina"
        String machineType = "NovaSeq"
        String center = "NYGC"

    }
    String libraryName = select_first([fastqsAlign.limsLibraryName, fastqsAlign.clientSampleId])
    command {
        minimap2 \
        -a \
        -xsr \
        -Y \
        -t ~{minimapThreads} \
        -R '@RG\tID:~{fastqsAlign.readGroupId}\tPL:~{platform}\tPM:~{machineType}\tLB:~{libraryName}\tDS:hg38\tSM:~{fastqsAlign.clientSampleId}\tCN:~{center}\tPU:${fastqsAlign.readGroupPlatformUnit}' \
        ~{bwaReference.fasta} \
        ~{fastqsAlign.fastqR1} \
        ~{fastqsAlign.fastqR2} \
        | samtools view \
        -@ ~{threads} \
        -Shb \
        -o ~{laneBamPath} \
        -
    }

    output {
        File laneBam = laneBamPath
    }

    runtime {
        mem: memoryGb + "G"
        cpus: totalThreads
        cpu : totalThreads
        memory : memoryGb + " GB"
        disks: "local-disk " + diskSize + " LOCAL"
        docker : "gcr.io/nygc-compbio/minimap2:2.20-r1061"
        runtime_minutes: "500"
    }
}

task ShortAlignMark {
    input {
        # command
        File laneBam
        String bamBase
        String laneBamMarkPath = "~{bamBase}.readgroup_mark.bam"
        # resources
        Int memoryGb = 16
        Int diskSize
        Int preemptible = 3
        Int timeLimit = 500
    }

    command {
        set -e -o pipefail
        timeout ~{timeLimit}m \
        filter_bam \
        -I ~{laneBam} \
        -A1 30 \
        -A2 30 \
        -o ~{laneBamMarkPath} \
        | samtools view \
        -b \
        -o ~{laneBamMarkPath} \
        -
    }

    output {
        File laneBamMark = laneBamMarkPath
    }

    runtime {
        mem: memoryGb + "G"
        memory : memoryGb + " GB"
        docker : "gcr.io/nygc-public/nygc-short-alignment-marking@sha256:8adc1ac0417d080cb6d6dd8901a5282c3cca497cff1b9800c2f51b0058b15bde"
        disks: "local-disk " + diskSize + " HDD"
        preemptible: preemptible
        runtime_minutes: "500"
    }
}

task Fixmate {
    input {
        #command
        File laneBamMark
        String bamBase
        String laneFixmateBamPath = "~{bamBase}.readgroup_fixmate.bam"
        # resources
        Int memoryGb = 8
        Int diskSize
        Int timeLimit = 500
    }

    Int jvmHeap = memoryGb * 750  # Heap size in Megabytes. mem is in GB. (75% of mem)

    command {
        timeout ~{timeLimit}m \
        gatk \
        FixMateInformation \
        --java-options "-Xmx~{jvmHeap}m -XX:ParallelGCThreads=1" \
        --MAX_RECORDS_IN_RAM 2000000 \
        --VALIDATION_STRINGENCY SILENT \
        --ADD_MATE_CIGAR true \
        --ASSUME_SORTED true \
        -I ~{laneBamMark} \
        -O ~{laneFixmateBamPath}
    }

    output {
        File laneFixmateBam = laneFixmateBamPath
    }

    runtime {
        mem: memoryGb + "G"
        memory : memoryGb + " GB"
        docker : "gcr.io/nygc-public/broadinstitute/gatk4@sha256:b3bde7bc74ab00ddce342bd511a9797007aaf3d22b9cfd7b52f416c893c3774c"
        disks: "local-disk " + diskSize + " HDD"
        runtime_minutes: "500"
    }
}
