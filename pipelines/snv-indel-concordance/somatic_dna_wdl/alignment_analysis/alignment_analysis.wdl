version 1.0

import "../wdl_structs.wdl"

task UnclippedPercentTrim {
    input {
        Bam finalBam
        IndexedReference referenceFa
        String sampleId
        File adaptersFa
        File teloTargetIndexGem # telo_target_index.gem
        
        String gemLengthPath = "~{sampleId}_total_telo_length.txt"
        Float maxMismatches = 0.2
        String alignmentHistoPath = "~{sampleId}.alignment.pdf"
        String alignmentCsvPath = "~{sampleId}.alignment.csv"
        String r1MappedFastqPath = "~{sampleId}.R1_mapped.fastq"
        String r2MappedFastqPath = "~{sampleId}.R2_mapped.fastq"
        String singMappedFastqPath = "~{sampleId}.singleton_mapped.fastq"
        String dollarSign = "$"
        File? serviceAccountKey
        String? gcpProject
        
        Int threads = 8
        Int samtoolsThreads = 2
        Int memoryGb = 10
        Int diskSize
    }
    
    command <<<
        set -e -o pipefail
        
        function gcs_auth {
            serviceAccountKey=~{serviceAccountKey}
            if [ -f "~{dollarSign}serviceAccountKey" ]; then
                export GOOGLE_APPLICATION_CREDENTIALS=~{serviceAccountKey}
                export GOOGLE_CLOUD_PROJECT=~{gcpProject}
                # expires in 60 min
                export GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token`
                export GCS_REQUESTER_PAYS_PROJECT=~{gcpProject}
            fi
        }
        
        # feed in reads skipping secondary or supplementary reads
        # and reads marked as vendor qc failed
        command="gcs_auth; samtools view \
        --reference ~{referenceFa.fasta} \
        -h \
        -F 2816 \
        --threads ~{threads} \
        ~{finalBam.bam} \
        | /pre_filter.py \
        | samtools fastq \
        --threads ~{samtoolsThreads} \
        - \
        | /trim_adapters.py \
        ~{adaptersFa} \
        ~{samtoolsThreads} \
        | gem-mapper \
        -T ~{threads} \
        --verbose \
        -I ~{teloTargetIndexGem} \
        -m ~{maxMismatches} \
        -e ~{maxMismatches} \
        --mismatch-alphabet ATCGN \
        --fast-mapping \
        -q ignore \
        | /describe_alignments.py \
        ~{alignmentHistoPath} \
        ~{alignmentCsvPath} \
        | /gem_to_fastq.py \
        ~{r1MappedFastqPath} \
        ~{r2MappedFastqPath} \
        ~{singMappedFastqPath} \
        ~{gemLengthPath}"
        
        eval ~{dollarSign}{command} || eval ~{dollarSign}{command}
    >>>
    
    output {
        File gemLength = "~{gemLengthPath}"
        File r2MappedFastq = "~{r2MappedFastqPath}"
        File r1MappedFastq = "~{r1MappedFastqPath}"
        File singMappedFastq =  "~{singMappedFastqPath}"
        File alignmentHisto = "~{alignmentHistoPath}"
        File gemLength = "~{gemLengthPath}"
        File alignmentCsv = "~{alignmentCsvPath}"

    }
    
    runtime {
        mem: memoryGb + "G"
        cpus: threads
        cpu : threads
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-comp-s-fd4e/telomeasures@sha256:6bd170351968ad3fe495c81202520258355f28ca7b7a4d3d77865ac9398c4b7c"
    }
    
    parameter_meta {
        finalBam: {
            description: "Input BAM filename",
            category: "required",
            localization_optional: true
        }
    }
}

task UnclippedPercent {
    input {
        Bam finalBam
        IndexedReference referenceFa
        String sampleId
        File teloTargetIndexGem # telo_target_index.gem
        
        String gemLengthPath = "~{sampleId}_total_telo_length.txt"
        Float maxMismatches = 0.2
        String alignmentHistoPath = "~{sampleId}.alignment.pdf"
        String alignmentCsvPath = "~{sampleId}.alignment.csv"
        String r1MappedFastqPath = "~{sampleId}.R1_mapped.fastq"
        String r2MappedFastqPath = "~{sampleId}.R2_mapped.fastq"
        String singMappedFastqPath = "~{sampleId}.singleton_mapped.fastq"
        String dollarSign = "$"
        File? serviceAccountKey
        String? gcpProject
        
        Int threads = 8
        Int samtoolsThreads = 2
        Int memoryGb = 10
        Int diskSize
    }
    
    command <<<
        set -e -o pipefail
        
        function gcs_auth {
            serviceAccountKey=~{serviceAccountKey}
            if [ -f "~{dollarSign}serviceAccountKey" ]; then
                export GOOGLE_APPLICATION_CREDENTIALS=~{serviceAccountKey}
                export GOOGLE_CLOUD_PROJECT=~{gcpProject}
                # expires in 60 min
                export GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token`
                export GCS_REQUESTER_PAYS_PROJECT=~{gcpProject}
            fi
        }
        
        # feed in reads skipping secondary or supplementary reads
        # and reads marked as vendor qc failed
        command="gcs_auth; samtools view \
        --reference ~{referenceFa.fasta} \
        -h \
        -F 2816 \
        --threads ~{threads} \
        ~{finalBam.bam} \
        | /pre_filter.py \
        | samtools fastq \
        --threads ~{samtoolsThreads} \
        - \
        | gem-mapper \
        -T ~{threads} \
        --verbose \
        -I ~{teloTargetIndexGem} \
        -m ~{maxMismatches} \
        -e ~{maxMismatches} \
        --mismatch-alphabet ATCGN \
        --fast-mapping \
        -q ignore \
        | /describe_alignments.py \
        ~{alignmentHistoPath} \
        ~{alignmentCsvPath} \
        | /gem_to_fastq.py \
        ~{r1MappedFastqPath} \
        ~{r2MappedFastqPath} \
        ~{singMappedFastqPath} \
        ~{gemLengthPath}"
        
        eval ~{dollarSign}{command} || eval ~{dollarSign}{command}
    >>>
    
    output {
        File gemLength = "~{gemLengthPath}"
        File r2MappedFastq = "~{r2MappedFastqPath}"
        File r1MappedFastq = "~{r1MappedFastqPath}"
        File singMappedFastq =  "~{singMappedFastqPath}"
        File alignmentHisto = "~{alignmentHistoPath}"
        File gemLength = "~{gemLengthPath}"
        File alignmentCsv = "~{alignmentCsvPath}"

    }
    
    runtime {
        mem: memoryGb + "G"
        cpus: threads
        cpu : threads
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-comp-s-fd4e/telomeasures@sha256:6bd170351968ad3fe495c81202520258355f28ca7b7a4d3d77865ac9398c4b7c"
    }
    
    parameter_meta {
        finalBam: {
            description: "Input BAM filename",
            category: "required",
            localization_optional: true
        }
    }
}

task TeloRatio {
    input {
        String sampleId
        File alignmentCsv 
        String telomeasuresPath = "~{sampleId}.telomeasures.summary.csv"
        
        Int threads = 1
        Int memoryGb = 1
        Int diskSize = 4
    }
    
    command {
        /telomeasures.py \
        --sample-id ~{sampleId} \
        --alignment-lengths ~{alignmentCsv} \
        > ~{telomeasuresPath}
    }
    
    output {
        File telomeasures = "~{telomeasuresPath}"
    }
    
    runtime {
        mem: memoryGb + "G"
        cpus: threads
        cpu : threads
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-comp-s-fd4e/telomeasures@sha256:6bd170351968ad3fe495c81202520258355f28ca7b7a4d3d77865ac9398c4b7c"
        
    }
}


task TelFusDetectorCaller {
    input {
        Int threads = 4
        Int memoryGb = 10
        Int diskSize = 10
        String sampleId
        String analysisId
        Bam finalBam
        String genome = "Hg38"
        String fusionsFilteredAutoPath = "~{sampleId}.summary_fusions.filtered.tsv"
        String fusionsPassAutoPath = "~{sampleId}.summary_fusions.pass.tsv"
        String allChromosomesCovAutoPath = "~{sampleId}.all_chromosomes.coverage.tsv"
        
        String fusionsFilteredPath = "~{analysisId}.summary_fusions.filtered.tsv"
        String fusionsPassPath = "~{analysisId}.summary_fusions.pass.tsv"
        String allChromosomesCovPath = "~{analysisId}.all_chromosomes.coverage.tsv"
        # needed for work on cloud without localizing
        File? serviceAccountKey
        String? gcpProject
        String dollarSign = "$"
    }
    command {
        gcpProject=~{gcpProject}
        serviceAccountKey=~{serviceAccountKey}
        if [  "~{dollarSign}gcpProject" != "" ]; then
            if [ -f "~{dollarSign}serviceAccountKey" ]; then
                export GOOGLE_APPLICATION_CREDENTIALS=~{serviceAccountKey}
            fi
            export GOOGLE_CLOUD_PROJECT=~{gcpProject}
            # expires in 60 min
            export GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token`
            export GCS_REQUESTER_PAYS_PROJECT=~{gcpProject}
        fi
        mkdir tmp/
        
        python3 \
        /TelFusDetectorCaller.py \
        --bam  ~{finalBam.bam} \
        --genome ~{genome} \
        --outfolder . \
        --sample ~{sampleId} \
        --tmpfolder tmp/ \
        --threads ~{threads}
        cp ~{fusionsFilteredAutoPath} \
        ~{fusionsFilteredAutoPath}.cp
        cp ~{fusionsFilteredAutoPath}.cp  \
        ~{fusionsFilteredPath}
        
        cp ~{fusionsPassAutoPath} \
        ~{fusionsPassAutoPath}.cp
        cp ~{fusionsPassAutoPath}.cp  \
        ~{fusionsPassPath}
        
        cp ~{allChromosomesCovAutoPath} \
        ~{allChromosomesCovAutoPath}.cp
        cp ~{allChromosomesCovAutoPath}.cp \
        ~{allChromosomesCovPath}
    }
    output {
        File fusionsFiltered = "~{fusionsFilteredPath}"
        File fusionsPass = "~{fusionsPassPath}"
        File allChromosomesCov = "~{allChromosomesCovPath}"
    }
    runtime {
        cpus: threads
        cpu : threads
        mem: memoryGb + "G"
        memory : memoryGb + "GB"
        disks: "local-disk " + diskSize + " HDD"
        docker : "gcr.io/nygc-comp-s-fd4e/telfusdetector@sha256:e7745f42c90164a19d47810620521f99cb2fe98628ceda2e284bd0bdbb1f2b83"
        runtime_minutes: "1440"
    }
    parameter_meta {
        finalBam: {
            description: "Input BAM filename",
            category: "required",
            localization_optional: true
        }
    }
}

task TelFusDetectorRates {
    input {
        Int memoryGb = 2
        Int diskSize = 2
        File fusionsPass
        Int purity = 1
        String sampleId
        String analysisId
        String genome = "Hg38"
        String fusionsRatesPath = "~{analysisId}.fusion_rates.tsv"
        String fusionsRatesAutoPath = "~{sampleId}.fusion_rates.tsv"
        Array[String] groupings = ["Chromosomal_region", "Orientation", "Filter"]
    }
    command {

        python3 \
        /TelFusDetectorRates.py \
        --fusion_file ~{fusionsPass} \
        --purity ~{purity} \
        --variables ~{sep=" " groupings} \
        --outfile ~{fusionsRatesPath}
        
        cp ~{fusionsRatesAutoPath} \
        ~{fusionsRatesAutoPath}.cp
        cp ~{fusionsRatesAutoPath}.cp \
        ~{fusionsRatesPath}
    }
    output {
        File fusionsRates = "~{fusionsRatesPath}"
    }
    runtime {
        mem: memoryGb + "G"
        memory : memoryGb + "GB"
        disks: "local-disk " + diskSize + " HDD"
        docker : "gcr.io/nygc-comp-s-fd4e/telfusdetector@sha256:e7745f42c90164a19d47810620521f99cb2fe98628ceda2e284bd0bdbb1f2b83"
        runtime_minutes: "120"
    }
}

task GetSampleNameFile {
    input {
        File finalBam
        File finalBai
        IndexedReference referenceFa
        String fileId
        String sampleIdPath = "sampleId.~{fileId}.txt"
        Int memoryGb = 1
        Int diskSize
    }

    Int jvmHeap = memoryGb * 750  # Heap size in Megabytes. mem is in GB. (75% of mem)
    command {
            /gatk/gatk \
            --java-options "-Xmx~{jvmHeap}m -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" \
            GetSampleName \
            -R ~{referenceFa.fasta} \
            -I ~{finalBam} \
            --read-index ~{finalBai} \
            -O ~{sampleIdPath}
    }

    output {
        File bamSampleId = "~{sampleIdPath}"
    }

    runtime {
        mem: memoryGb + "G"
        memory : memoryGb + "GB"
        disks: "local-disk " + diskSize + " HDD"
        docker: "gcr.io/nygc-public/broadinstitute/gatk4@sha256:b3bde7bc74ab00ddce342bd511a9797007aaf3d22b9cfd7b52f416c893c3774c"
        runtime_minutes: "60"
    }
    
    parameter_meta {
        finalBam: {
            localization_optional: true
        }
        
        finalBai: {
            localization_optional: true
        }
        
        referenceFa: {
            localization_optional: true
        }
    }
}

task GetSampleName {
    input {
        File finalBam
        File finalBai
        IndexedReference referenceFa
        String sampleIdPath = "sampleId.txt"
        Int memoryGb = 1
        Int diskSize
    }

    Int jvmHeap = memoryGb * 750  # Heap size in Megabytes. mem is in GB. (75% of mem)
    command {
            /gatk/gatk \
            --java-options "-Xmx~{jvmHeap}m -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" \
            GetSampleName \
            -R ~{referenceFa.fasta} \
            -I ~{finalBam} \
            --read-index ~{finalBai} \
            -O ~{sampleIdPath}
    }

    output {
        String bamSampleId = read_string("~{sampleIdPath}")
    }

    runtime {
        mem: memoryGb + "G"
        memory : memoryGb + "GB"
        disks: "local-disk " + diskSize + " HDD"
        docker: "gcr.io/nygc-public/broadinstitute/gatk4@sha256:b3bde7bc74ab00ddce342bd511a9797007aaf3d22b9cfd7b52f416c893c3774c"
        runtime_minutes: "60"
    }
    
    parameter_meta {
        finalBam: {
            localization_optional: true
        }
        
        finalBai: {
            localization_optional: true
        }
        
        referenceFa: {
            localization_optional: true
        }
    }
}

task UpdateBamSampleName {
    input {
        Bam finalBam
        IndexedReference referenceFa
        String sampleId
        String outputPrefix = "~{sampleId}"
        String headerPath = "~{outputPrefix}.reheader.txt"
        String reheaderBamPath = "~{outputPrefix}.reheader.bam"
        String bamIndexPath = "~{outputPrefix}.reheader.bai"
        # resources
        Int diskSize
        Int threads = 4
        Int memoryGb = 8
    }

    command {
    
        set -e -o pipefail
        
        samtools \
        view -H \
        ~{finalBam.bam} \
        | sed "/^@RG/ s/SM:\S*/SM:~{sampleId}/" \
        > ~{headerPath}
        
        # also convert to bam (incase it is a cram file)
        samtools \
        reheader \
        ~{headerPath} \
        ~{finalBam.bam} \
        | samtools \
        view \
        -b \
        -T ~{referenceFa.fasta} \
        -t ~{referenceFa.index} \
        -o ~{reheaderBamPath} \
        --verbosity=8 \
        --threads ~{threads} \
        -
        
        samtools \
        index \
        -@ ~{threads} \
        ~{reheaderBamPath} \
        ~{bamIndexPath}
    }

    output {
        Bam reheaderBam = object {
            bam : reheaderBamPath,
            bamIndex : bamIndexPath
        }
    }

    runtime {
        docker : "gcr.io/nygc-public/samtools@sha256:32f29fcd7af01b3941e6f93095e8d899741e81b50bcc838329bd8df43e120cc3"
        disks: "local-disk " + diskSize + " LOCAL"
        memory: memoryGb + "GB"
        mem: memoryGb + "G"
        cpu: threads
        cpus: threads
        runtime_minutes: "90"
    }
}



task BedtoolsIntersect {
    input {
        String mantisBedByIntervalListPath
        File mantisBed
        File intervalListBed
        Int memoryGb = 1
    }

    command {
        bedtools \
        intersect \
        -a ~{mantisBed} \
        -b ~{intervalListBed} \
        -u \
        > ~{mantisBedByIntervalListPath} \
    }

    output {
        File mantisBedByIntervalList = "~{mantisBedByIntervalListPath}"
    }

    runtime {
        mem: memoryGb + "G"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/bedtools@sha256:9e737f5c96c00cf3b813d419d7a7b474c4013c9aa9dfe704eb36417570c6474e"
        runtime_minutes: "60"
    }
}


task MantisExome {
    input {
        String pairName
        String mantisExomeTxtPath = "~{pairName}.mantis.WGS-targeted.txt"
        String mantisWxsKmerCountsPath = "~{pairName}.mantis.WGS-targeted.kmer_counts.txt"
        
        String mantisWxsKmerCountsFilteredPath = "~{pairName}.mantis.WGS-targeted.kmer_counts_filtered.txt"
        String mantisWxsStatusPath = "~{pairName}.mantis.WGS-targeted.txt.status"
        
        Bam tumorFinalBam
        Bam normalFinalBam
        String tumorFinalBamPath = basename(tumorFinalBam.bam)
        String tumorFinalBamIndexPath = basename(tumorFinalBam.bamIndex)
        String normalFinalBamPath = basename(normalFinalBam.bam)
        String normalFinalBamIndexPath = basename(normalFinalBam.bamIndex)

        File mantisBedByIntervalList
        IndexedReference referenceFa
        Int threads = 8
        Int memoryGb = 4
        Int diskSize = ceil( size(tumorFinalBam.bam, "GB") + size(normalFinalBam.bam, "GB")) + 30

    }

    command {
        set -e -o pipefail

        # make a .bam.bai index available
        # normal
        ln -s \
        ~{normalFinalBam.bam} \
        ~{normalFinalBamPath}

        ln -s \
        ~{normalFinalBam.bamIndex} \
        ~{normalFinalBamIndexPath}

        ln -s \
        ~{normalFinalBamIndexPath} \
        ~{normalFinalBamPath}.bai

        # tumor
        ln -s \
        ~{tumorFinalBam.bam} \
        ~{tumorFinalBamPath}

        ln -s \
        ~{tumorFinalBam.bamIndex} \
        ~{tumorFinalBamIndexPath}

        ln -s \
        ~{tumorFinalBamIndexPath} \
        ~{tumorFinalBamPath}.bai

        ls -thl

        python \
        /MANTIS-1.0.4/mantis.py \
        --bedfile ~{mantisBedByIntervalList} \
        --genome ~{referenceFa.fasta} \
        -mrq 20.0 \
        -mlq 25.0 \
        -mlc 20 \
        -mrr 1 \
        --threads ~{threads} \
        -n ~{normalFinalBamPath} \
        -t ~{tumorFinalBamPath} \
        -o ~{mantisExomeTxtPath}
        
        # create blank mock file if
        # no status was called because
        # regions failed mantis QC thresholds
        
        if [[ -e ~{mantisWxsStatusPath} ]]; then
            echo 'Mantis status called'
        else
            echo 'Mantis status not called (too many regions filtered out)'
            echo '' > ~{mantisWxsStatusPath}
        fi

    }

    output {
        File mantisWxsKmerCountsFinal = "~{mantisWxsKmerCountsPath}"
        File mantisWxsKmerCountsFiltered = "~{mantisWxsKmerCountsFilteredPath}"
        File mantisWxsStatus = "~{mantisWxsStatusPath}"
        File mantisExomeTxt = "~{mantisExomeTxtPath}"
    }

    runtime {
        mem: memoryGb + "G"
        cpus: threads
        cpu : threads
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/mantis@sha256:9cf1311c5198b8fa5fecff387a50dfa9408707f7b914a99dc548c6eb14f42c19"
        runtime_minutes: "500"
    }
}

task MantisRethreshold {
    input {
        String pairName
        String mantisStatusFinalPath = "~{pairName}.mantis.WGS-targeted.status.final.tsv"
        String normal
        File mantisWxsStatus
    }

    command {
        set -e -o pipefail
        
        if [[ -s ~{mantisWxsStatus} ]]; then
            python \
            /reset_mantis.py \
            ~{mantisWxsStatus} \
            ~{mantisStatusFinalPath} \
            ~{normal}
        else
            echo \
            'MSI_Status    MSI_Score    Threshold' \
            > ~{mantisStatusFinalPath}
        fi

    }

    output {
        File mantisStatusFinal = "~{mantisStatusFinalPath}"
    }

    runtime {
        docker : "gcr.io/nygc-public/somatic_dna_tools@sha256:f281e73ddf515f3a5db6e766ef12fb2331dafa2b27c88e426764dcd8b4a7e19b"
        runtime_minutes: "60"
    }
}

task GetChr6Contigs {
    input  {
        IndexedReference referenceFa
        Int diskSize
        Int memoryGb = 2
    }

    command {
        /lookup_contigs.py ~{referenceFa.fasta}
    }

    output {
        String chr6Contigs = read_string(stdout())
    }

    runtime {
        mem: memoryGb + "G"
        memory : memoryGb + "GB"
        disks: "local-disk " + diskSize + " HDD"
        docker : "gcr.io/nygc-public/hla_prep@sha256:a490cf449eeb98b997f0dd87ff1c23ff77d724c7c2072b6c44f75a713ecc2d36"
        runtime_minutes: "500"
    }
}

task GemSelect {
    input {
        Int threads = 8
        Int samtoolsThreads = 8
        Int gemThreads = 8
        Int memoryGb = 4
        Int diskSize
        String sampleId
        String chr6Contigs
        Bam finalBam
        File kouramiFastaGem1Index
        String r1FilePath = "~{sampleId}.first_pair"
        String r2FilePath = "~{sampleId}.second_pair"
        Float maxMismatches = 0.04
        String alignmentHistoPath = "~{sampleId}.alignment.pdf"
        String r1MappedFastqPath = "~{sampleId}.R1_mapped.fastq"
        String r2MappedFastqPath = "~{sampleId}.R2_mapped.fastq"
    }

    command {
        set -e -o pipefail

        samtools view \
        --threads ~{samtoolsThreads} \
        -h \
        -f 1 \
        ~{finalBam.bam} \
        ~{chr6Contigs} \
        | /note_pair.py \
        ~{r1FilePath} \
        ~{r2FilePath} \
        | samtools fastq \
        --threads ~{samtoolsThreads} \
        - \
        | gem-mapper \
        -T ~{gemThreads} \
        --verbose \
        -I ~{kouramiFastaGem1Index} \
        -m ~{maxMismatches} \
        -e ~{maxMismatches} \
        --mismatch-alphabet ATCGN \
        --fast-mapping \
        -q ignore \
        | /describe_alignments.py \
        ~{alignmentHistoPath} \
        | /gem_to_fastq.py \
        ~{r1MappedFastqPath} \
        ~{r2MappedFastqPath}
    }

    output {
        File r2File = "~{r1FilePath}"
        File r2MappedFastq = "~{r2MappedFastqPath}"
        File r1File = "~{r2FilePath}"
        File r1MappedFastq = "~{r1MappedFastqPath}"
        File alignmentHisto = "~{sampleId}.alignment.pdf"
    }

    runtime {
        mem: memoryGb + "G"
        cpus: threads
        cpu : threads
        memory : memoryGb + "GB"
        disks: "local-disk " + diskSize + " HDD"
        docker : "gcr.io/nygc-public/hla_prep@sha256:a490cf449eeb98b997f0dd87ff1c23ff77d724c7c2072b6c44f75a713ecc2d36"
        runtime_minutes: "500"
    }
}

task LookUpMates {
    input {
        Int memoryGb = 2
        Int diskSize = 4
        String sampleId
        String r1UnmappedFilePath = "~{sampleId}.first_pair_unmapped"
        String r2UnmappedFilePath = "~{sampleId}.second_pair_unmapped"
        File r2File
        File r2MappedFastq
        File r1File
        File r1MappedFastq

    }

    command {
        /look_up_mates.py \
        ~{r1File} \
        ~{r2File} \
        ~{r1MappedFastq} \
        ~{r2MappedFastq} \
        ~{r1UnmappedFilePath} \
        ~{r2UnmappedFilePath}
    }

    output {
        File r1UnmappedFile = "~{r1UnmappedFilePath}"
        File r2UnmappedFile = "~{r2UnmappedFilePath}"
    }

    runtime {
        mem: memoryGb + "G"
        memory : memoryGb + "GB"
        disks: "local-disk " + diskSize + " HDD"
        docker : "gcr.io/nygc-public/hla_prep@sha256:a490cf449eeb98b997f0dd87ff1c23ff77d724c7c2072b6c44f75a713ecc2d36"
        runtime_minutes: "60"
    }
}

task GetMatesLocalize {
    input {
        Int threads = 8
        Int samtoolsThreads = 4
        Int memoryGb = 2
        Int diskSize
        String sampleId
        String r1UnmappedFastqPath = "~{sampleId}.R1_unmapped.fastq"
        String r2UnmappedFastqPath = "~{sampleId}.R2_unmapped.fastq"
        Bam finalBam
        File r1UnmappedFile
        File r2UnmappedFile
    }

    command {
        set -e -o pipefail

        samtools view \
        --threads ~{samtoolsThreads} \
        -h \
        -f 1 \
        ~{finalBam.bam} \
        | python3 /get_mates.py \
        ~{r1UnmappedFile} \
        ~{r2UnmappedFile} \
        | samtools fastq \
        --threads ~{samtoolsThreads} \
        -1 ~{r1UnmappedFastqPath} \
        -2 ~{r2UnmappedFastqPath} \
        -
    }

    output {
        File r1UnmappedFastq = "~{r1UnmappedFastqPath}"
        File r2UnmappedFastq = "~{r2UnmappedFastqPath}"
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
    
    parameter_meta {
        finalBam: {
            description: "Input BAM filename",
            category: "required",
            localization_optional: false
        }
    }
}

task GetMates {
    input {
        Int threads = 8
        Int samtoolsThreads = 4
        Int memoryGb = 2
        Int diskSize
        String sampleId
        String r1UnmappedFastqPath = "~{sampleId}.R1_unmapped.fastq"
        String r2UnmappedFastqPath = "~{sampleId}.R2_unmapped.fastq"
        Bam finalBam
        File r1UnmappedFile
        File r2UnmappedFile
        
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
        --threads ~{samtoolsThreads} \
        -h \
        -f 1 \
        ~{finalBam.bam} \
        | python3 /get_mates.py \
        ~{r1UnmappedFile} \
        ~{r2UnmappedFile} \
        | samtools fastq \
        --threads ~{samtoolsThreads} \
        -1 ~{r1UnmappedFastqPath} \
        -2 ~{r2UnmappedFastqPath} \
        -
    }

    output {
        File r1UnmappedFastq = "~{r1UnmappedFastqPath}"
        File r2UnmappedFastq = "~{r2UnmappedFastqPath}"
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
    
    parameter_meta {
        finalBam: {
            description: "Input BAM filename",
            category: "required",
            localization_optional: true
        }
    }
}

task SortFastqs {
    input {
        Int memoryGb = 2
        Int diskSize = 4
        String sampleId
        String fastqPairId
        String sortedFastqPath = "~{sampleId}.~{fastqPairId}_sorted.fastq"
        File chr6MappedFastq
        File chr6MappedMatesFastq
    }

    command {
        set -e -o pipefail

        cat \
        ~{chr6MappedFastq} \
        ~{chr6MappedMatesFastq} \
        | seqkit fx2tab \
        | /match_header.py \
        | sort \
        --dictionary-order \
        -k1,1 \
        -T \
        . \
        | seqkit tab2fx \
        > ~{sortedFastqPath}
    }

    output {
        File sortedFastq = "~{sortedFastqPath}"
    }

    runtime {
        mem: memoryGb + "G"
        memory : memoryGb + "GB"
        disks: "local-disk " + diskSize + " HDD"
        docker : "gcr.io/nygc-public/hla_prep@sha256:a490cf449eeb98b997f0dd87ff1c23ff77d724c7c2072b6c44f75a713ecc2d36"
        runtime_minutes: "60"
    }
}

task AlignToPanel {
    input {
        Int threads = 8
        Int bwaThreads = 4
        Int samtoolsSortThreads = 4
        Int memoryGb = 4
        Int diskSize = 4
        String sampleId
        String kouramiBamPath = "~{sampleId}.kourami.bam"
        File r2SortedFastq
        # mergedHlaPanel
        BwaReference kouramiReference
        File r1SortedFastq
    }

    command {
        set -e -o pipefail

        bwa mem \
        -t ~{bwaThreads} \
        ~{kouramiReference.fasta} \
        ~{r1SortedFastq} \
        ~{r2SortedFastq} \
        | samtools sort \
        --threads ~{samtoolsSortThreads} \
        -m 10G \
        -o ~{kouramiBamPath}
    }

    output {
        File kouramiBam = "~{kouramiBamPath}"
    }

    runtime {
        mem: memoryGb + "G"
        cpus: threads
        cpu : threads
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/bwa-kit@sha256:0642151a32fe8f90ece70cde3bd61a03c7421314c37c1de2c0ee5e368d2bfc7a"
        runtime_minutes: "60"
    }
}

task Kourami {
    input {
        Int threads = 1
        Int memoryGb = 8
        String sampleId
        File kouramiBam
        String resultPrefix = "~{sampleId}.kourami"
    }

    Int jvmHeap = memoryGb * 750  # Heap size in Megabytes. mem is in GB. (75% of mem)

    command {
        set -e -o pipefail
        
        java \
        -Xmx~{jvmHeap}m -XX:ParallelGCThreads=4 \
        -jar /Kourami.jar \
        -d /kourami-0.9.6/db/ \
        -o ~{resultPrefix} \
        ~{kouramiBam}
        
    }

    output {
        File result = "~{resultPrefix}.result"
    }

    runtime {
        mem: memoryGb + "G"
        cpus: threads
        cpu : threads
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/kourami@sha256:d4b906b979c24ee4669fdbf7ee1dfbdeb5c89d0e34b4b4aaf21ee070e988d74b"
        runtime_minutes: "60"
    }
}

task Angsd {
    input {
        Bam normalFinalBam
        File fastNgsAdmixSites
        File fastNgsAdmixSitesBin
        File fastNgsAdmixSitesIdx
        File fastNgsAdmixChroms
        Int threads
        String outprefix
        Int memoryGb = 25
        Int diskSize = ceil(size(normalFinalBam.bam, "GB")) + 30
    }

    command {
        set -exo pipefail

        angsd \
            -i ~{normalFinalBam.bam} \
            -GL 2 \
            -rf ~{fastNgsAdmixChroms} \
            -sites ~{fastNgsAdmixSites} \
            -doMajorMinor 3 \
            -doGlf 2 \
            -minMapQ 30 \
            -minQ 20 \
            -doDepth 1 \
            -doCounts 1 \
            -nThreads ~{threads} \
            -out ~{outprefix}
    }

    runtime {
        memory: memoryGb + "G"
        mem: memoryGb + "G"
        cpu: threads
        cpus: threads
        docker: "gcr.io/nygc-public/angsd@sha256:cd13820de0bc8d400c3e3ff96be6b885b6d3289d53fda56de30bd08508a0bac7"
        disks: "local-disk " + diskSize + " HDD"
        runtime_minutes: "500"
    }

    output {
        File beagleFile = "${outprefix}.beagle.gz"
        File beagleLog = "${outprefix}.arg"
        File beagleDepth = "${outprefix}.depthGlobal"
        File beagleSample = "${outprefix}.depthSample"
    }
}

task FastNgsAdmix {
    input {
        File beagleFile
        File fastNgsAdmixRef
        File fastNgsAdmixNind
        String outprefix
        Int memoryGb = 15
        Int threads = 1
        Int diskSize = 30
    }

    command {
        set -exo pipefail

        fastNGSadmix  \
            -likes ~{beagleFile} \
            -fname ~{fastNgsAdmixRef} \
            -Nname ~{fastNgsAdmixNind} \
            -whichPops all \
            -out  ~{outprefix}
    }

    runtime {
        memory: memoryGb + "G"
        mem: memoryGb + "G"
        cpu: threads
        cpus: threads
        docker: "gcr.io/nygc-public/fastngsadmix@sha256:f0a336e9f193ab1b4f1484cbec56e1abef063913102d3126f4b3a6ed7784d7f1"
        disks: "local-disk " + diskSize + " HDD"
        runtime_minutes: "60"
    }

    output {
        File fastNgsAdmixQopt = "${outprefix}.qopt"
        File fastNgsAdmixLog = "${outprefix}.log"
    }
}
