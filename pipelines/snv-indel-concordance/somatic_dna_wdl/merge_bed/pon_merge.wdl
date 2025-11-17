version 1.0

import "../wdl_structs.wdl"

task MakeSinglePon {
    input {
        Int threads = 1
        Int memoryGb = 40
        Int diskSize
        String sampleId
        File vcfAnnotatedVep
        String simpleTable = "~{sampleId}.pon.table.bed"
        String singlePonPath = "~{sampleId}.pon.bed"
        String openBraket = "{"
        String closeBraket = "}"
    }
    
    command {
        echo "Start header..."
        cat ~{vcfAnnotatedVep} \
        | awk '~{openBraket}if(/^#/)print;else exit~{closeBraket}' \
        | grep "^#CHROM" \
        | cut -f 1-7 \
        > ~{simpleTable}
        
        set -e -o pipefail
        
        echo "Finish table..."
        grep -v "^#" ~{vcfAnnotatedVep} \
        | cut -f 1-7 \
        >> ~{simpleTable}
        
        echo "Refine table..."
        python \
        /count_pon_sites.py \
        ~{simpleTable} \
        ~{singlePonPath}
        
    }

    output {
        File singlePon = "~{singlePonPath}"
    }

    runtime {
        mem: memoryGb + "G"
        cpus: threads
        cpu : threads
        memory : memoryGb + "GB"
        disks: "local-disk " + diskSize + " HDD"
        docker : "gcr.io/nygc-public/somatic_dna_tools@sha256:f281e73ddf515f3a5db6e766ef12fb2331dafa2b27c88e426764dcd8b4a7e19b"
        runtime_minutes: "500"
    }
}

task MergePons {
    input {
        Int threads = 1
        Int memoryGb = 40
        Int diskSize = 300
        String name = "WGS_1000g_GRCh38"
        Array[File] singlePons
        String ponPath = "~{name}.pon.bed"
    }
    
    command {        
        python \
        /merge_pon_sites.py \
        ~{ponPath} \
        ~{sep=" " singlePons}
        
    }

    output {
        File pon = "~{ponPath}"
    }

    runtime {
        mem: memoryGb + "G"
        cpus: threads
        cpu : threads
        memory : memoryGb + "GB"
        disks: "local-disk " + diskSize + " HDD"
        docker : "gcr.io/nygc-public/somatic_dna_tools@sha256:f281e73ddf515f3a5db6e766ef12fb2331dafa2b27c88e426764dcd8b4a7e19b"
        runtime_minutes: "500"
    }
}

