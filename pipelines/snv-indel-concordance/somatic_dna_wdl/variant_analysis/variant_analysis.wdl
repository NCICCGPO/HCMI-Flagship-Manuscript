version 1.0

import "../wdl_structs.wdl"


task SigProfilerMatrixGenerator {
    input {
        String pairId
        String outputPrefix = "~{pairId}.cosmic.v3.4"
        String outputSbsPath = "~{pairId}.cosmic.v3.4.SBS96.all"
        String outputDbsPath = "~{pairId}.cosmic.v3.4.DBS78.all"
        String outputIdPath = "~{pairId}.cosmic.v3.4.ID83.all"
        # from FilterHighConfidence
        File filteredVcf
        String vepGenomeBuild
        Int memoryGb = 8
        Int diskSize = 16
    }
    command {
        set -e -o pipefail
        python \
        /SigMatrixGen.py \
        --vcf ~{filteredVcf} \
        --genome-build ~{vepGenomeBuild} \
        --output-prefix ~{outputPrefix}
        # copy from home to working dir
        cp ~/*all .
        ls -thl 
    }

    output {
       File outputSbs = "~{outputSbsPath}"
       File outputDbs = "~{outputDbsPath}"
       File outputId = "~{outputIdPath}"
    }

    runtime {
        mem: memoryGb + "G"
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-comp-s-fd4e/gcs_htslib_suite@sha256:47eba58683641905a58f31c05605c953dc1d888466e8d434a6ceff47b76df03e"
        runtime_minutes: "500"
    }
}

task MusicalRefit {
    input {
        String pairId
        String outputPrefix = "~{pairId}.cosmic.v3.4"
        File outputSbs
        File outputDbs
        File outputId
        String sigsSbsPath = "~{pairId}.cosmic.v3.4.musical.SBS.txt"
        String sigsDbsPath = "~{pairId}.cosmic.v3.4.musical.DBS.txt"
        String sigsIdPath = "~{pairId}.cosmic.v3.4.musical.ID.txt"
        File cosmicSbsReference
        File cosmicDbsReference
        File cosmicIdReference
        # from FilterHighConfidence
        Int memoryGb = 8
        Int diskSize = 16
    }
    command {
        python \
        /run_musical_refit.py \
        --input-sbs ~{outputSbs} \
        --input-dbs ~{outputDbs} \
        --input-id ~{outputId} \
        --sbs ~{cosmicSbsReference} \
        --dbs ~{cosmicDbsReference} \
        --id ~{cosmicIdReference} \
        --output-prefix ~{outputPrefix}
    }

    output {
       File sigsSbs = "~{sigsSbsPath}"
       File? sigsDbs = "~{sigsDbsPath}"
       File? sigsId = "~{sigsIdPath}"
    }

    runtime {
        mem: memoryGb + "G"
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-comp-s-fd4e/gcs_htslib_suite@sha256:47eba58683641905a58f31c05605c953dc1d888466e8d434a6ceff47b76df03e"
        runtime_minutes: "500"
    }
}

task Deconstructsig {
    input {
        String pairId
        String outputPrefix = "~{pairId}.cosmic.v3.2.deconstructSigs.signatures.highconfidence"
        File mainVcf
        String vepGenomeBuild
        String highConf = "TRUE"
        File cosmicSigs
        Int memoryGb = 8
        Int diskSize = 1
    }

    command {
        Rscript \
        /run_deconstructSigs.R \
        --highconf ~{highConf} \
        --file ~{mainVcf} \
        --ref ~{vepGenomeBuild} \
        --cosmic ~{cosmicSigs} \
        --output ~{outputPrefix} \
        --samplename ~{pairId}
    }

    output {
       File sigs = "~{outputPrefix}.txt"
       File counts = "~{outputPrefix}.counts.txt"
       File sigInput = "~{outputPrefix}.input.txt"
       File reconstructed = "~{outputPrefix}.reconstructed.txt"
       File diff = "~{outputPrefix}.diff.txt"
    }

    runtime {
        mem: memoryGb + "G"
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/deconstructsigs@sha256:009ddb6ed3ec2a0290a88b1e7027dd3caac1a2f5f3df3e8f68f410481d9323a3"
        runtime_minutes: "60"
    }
}
