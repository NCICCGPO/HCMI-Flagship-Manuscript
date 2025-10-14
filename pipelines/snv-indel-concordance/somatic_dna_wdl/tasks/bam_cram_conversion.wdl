version 1.0

import "../wdl_structs.wdl"

task SamtoolsBamToCram {
    input {
        # command
        Bam inputBam
        IndexedReference referenceFa
        String sampleId
        String sampleAnalysisId
        String cramPath = "~{sampleAnalysisId}.final.cram"
        String indexPath = "~{sampleAnalysisId}.final.cram.crai"
        # resources
        Int diskSize
        Int threads = 8
        Int memoryGb = 8
    }

    command {
        set -e -o pipefail

        set -x

        samtools \
        view \
        -C \
        -T ~{referenceFa.fasta} \
        -t ~{referenceFa.index} \
        -o ~{cramPath} ~{inputBam.bam} \
        --verbosity=8 \
        --threads ~{threads}

        # Reheader the cram with the http reference
        samtools view -H ~{cramPath} | sed "s%UR:~{referenceFa.fasta}%UR:~{referenceFa.httpFasta}%" > new_header.txt
        samtools reheader -i --no-PG new_header.txt ~{cramPath}

        # Index cram
        samtools index ~{cramPath} ~{indexPath}
    }

    output {
        SampleCramInfo cramInfo = object {
            sampleId : "~{sampleId}",
            finalCram : {
                            'cram' : "~{cramPath}",
                            'cramIndex' : "~{indexPath}"
                        }
             }
    }

    runtime {
        docker : "gcr.io/nygc-comp-s-fd4e/gcs_htslib_suite@sha256:47eba58683641905a58f31c05605c953dc1d888466e8d434a6ceff47b76df03e"
        disks: "local-disk " + diskSize + " LOCAL"
        memory: memoryGb + "GB"
        mem: memoryGb + "G"
        cpu: threads
        cpus: threads
        max_retries: 3
        runtime_minutes: "500"
    }
}

# note: there are possible issues converting cram to bam directly. May need cram->sam->bam
# https://github.com/gatk-workflows/seq-format-conversion/blob/master/cram-to-bam.wdl
# The reason this approach was chosen instead of converting CRAM to BAM directly using Samtools is because Samtools 1.3
# produces incorrect bins due to an old version of htslib included in the package.
# Samtools versions 1.4 & 1.5 have an NM issue that causes them to not validate with Picard.

task SamtoolsCramToBam {
    input {
        # command
        Cram inputCram
        IndexedReference referenceFa
        String sampleId
        String bamPath = "~{sampleId}.final.bam"
        String indexPath = "~{sampleId}.final.bai"
        # resources
        Int diskSize
        Int threads = 8
        Int memoryGb = 8
    }

    command {
        set -e -o pipefail

        samtools \
        view \
        -b \
        -T ~{referenceFa.fasta} \
        -t ~{referenceFa.index} \
        -o ~{bamPath} ~{inputCram.cram} \
        --verbosity=8 \
        --threads ~{threads}

        samtools index ~{bamPath} ~{indexPath}
    }

    output {
        SampleBamInfo bamInfo = object {
            sampleId: "~{sampleId}",
            finalBam: {
                           'bam' : "~{bamPath}",
                           'bamIndex' : "~{indexPath}"
                       }
            }
    }

    runtime {
        docker : "gcr.io/nygc-comp-s-fd4e/gcs_htslib_suite@sha256:47eba58683641905a58f31c05605c953dc1d888466e8d434a6ceff47b76df03e"
        disks: "local-disk " + diskSize + " HDD"
        memory: memoryGb + "GB"
        mem: memoryGb + "G"
        cpu: threads
        cpus: threads
        runtime_minutes: "500"
    }
}


task UniqueBams {
    input {
        File pairInfosJson
    }

    command {
        python /unique_bams.py \
        --pair-infos ~{pairInfosJson}
    }

    output {
        Array[SampleBamInfo] uniqueBams = read_json('unique_bams.json')
    }

    runtime {
        docker : "gcr.io/nygc-public/somatic_dna_tools@sha256:f281e73ddf515f3a5db6e766ef12fb2331dafa2b27c88e426764dcd8b4a7e19b"
        disks: "local-disk 10 HDD"
        memory: "1GB"
        mem: "1G"
        cpu: "1"
        cpus: "1"
        runtime_minutes: "60"
    }
}

task UpdateCramInfos {
    input {
        File pairInfosJson
        File normalInfosJson
        File cramInfosJson
    }

    command {
        python /update_cram_info.py \
        --pair-infos ~{pairInfosJson} \
        --normal-infos ~{normalInfosJson} \
        --cram-infos ~{cramInfosJson}
    }

    output {
        Array[SampleCramInfo] normalSampleCramInfos = read_json('normal_cram_infos.json')
        Array[PairCramInfo] pairCramInfos = read_json('pair_cram_infos.json')
    }

    runtime {
        docker : "gcr.io/nygc-public/somatic_dna_tools@sha256:f281e73ddf515f3a5db6e766ef12fb2331dafa2b27c88e426764dcd8b4a7e19b"
        disks: "local-disk 10 HDD"
        memory: "1GB"
        mem: "1G"
        cpu: "1"
        cpus: "1"
        runtime_minutes: "90"
    }

}
