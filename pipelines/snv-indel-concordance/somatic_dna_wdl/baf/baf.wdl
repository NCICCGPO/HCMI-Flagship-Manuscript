version 1.0

import "../wdl_structs.wdl"

task FilterForHetSnps {
    input {
        String sampleId
        String hetVcfPath = "~{sampleId}.haplotypecaller.gatk.het.vcf"
        String sellectionString = "'vc.getGenotype(\"~{sampleId}\").isHet()'"
        IndexedReference referenceFa
        # require file!
        # marked as optional so that pipeline can be dependent on input that may not
        # pass QC
        # do not run with out input file!
        File? germlineVcf

        Int memoryGb = 24
        Int diskSize = (ceil( size(germlineVcf, "GB") )  * 2 ) + 10
    }

    Int jvmHeap = memoryGb * 750  # Heap size in Megabytes. mem is in GB. (75% of mem)

    command {
        gatk \
        SelectVariants \
        --java-options "-Xmx~{jvmHeap}m" \
        -restrict-alleles-to BIALLELIC \
        -select-type SNP \
        -select ~{sellectionString} \
        -R ~{referenceFa.fasta} \
        --exclude-filtered \
        -V ~{germlineVcf} \
        -O ~{hetVcfPath}
    }

    output {
        File hetVcf = "~{hetVcfPath}"
    }

    runtime {
        mem: memoryGb + "G"
        memory : memoryGb + "GB"
        disks: "local-disk " + diskSize + " HDD"
        docker: "gcr.io/nygc-public/broadinstitute/gatk4@sha256:b3bde7bc74ab00ddce342bd511a9797007aaf3d22b9cfd7b52f416c893c3774c"
        runtime_minutes: "360"
    }
}

task FilterBaf {
    input {
        String sampleId
        String knownHetVcfPath = "~{sampleId}.haplotypecaller.gatk.known.het.vcf"
        File hetVcf
        Int memoryGb = 4
        Int diskSize = (ceil( size(hetVcf, "GB") )  * 2 ) + 10
    }

    command {
        python \
        /filter_baf.py \
        ~{hetVcf} \
        ~{knownHetVcfPath}
    }

    output {
        File knownHetVcf = "~{knownHetVcfPath}"
    }

    runtime {
        mem: memoryGb + "G"
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/somatic_dna_tools@sha256:f281e73ddf515f3a5db6e766ef12fb2331dafa2b27c88e426764dcd8b4a7e19b"
        runtime_minutes: "180"
    }
}

task AlleleCounts {
    input {
        String pairName
        String alleleCountsTxtPath = "~{pairName}.haplotypecaller.gatk.alleles.txt"
        IndexedReference referenceFa
        Bam normalFinalBam
        File knownHetVcf
        Bam tumorFinalBam
        Int memoryGb = 64

        Int diskSize = (ceil( size(knownHetVcf, "GB") )  * 2 ) + 10
        
        Int timeLimit = 600
        
        # needed for work on cloud without localizing
        File? serviceAccountKey
        String? gcpProject
    }

    command {
        serviceAccountKey=~{serviceAccountKey}
        if [ -f "$serviceAccountKey" ]; then
            export GOOGLE_APPLICATION_CREDENTIALS=~{serviceAccountKey}
            export GOOGLE_CLOUD_PROJECT=${gcpProject}
            # expires in 60 min
            export GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token`
            export GCS_REQUESTER_PAYS_PROJECT=~{gcpProject}
        fi
        
        timeout ~{timeLimit}m \
        python3 \
        /parse_bam_generate_features.py \
        --tumor_bam ~{tumorFinalBam.bam} \
        --normal_bam ~{normalFinalBam.bam} \
        --vcf ~{knownHetVcf} \
        --output ~{alleleCountsTxtPath} \
        --reference ~{referenceFa.fasta}
    }

    output {
        File alleleCountsTxt = "~{alleleCountsTxtPath}"
    }

    runtime {
        mem: memoryGb + "G"
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-comp-s-fd4e/gcs_htslib_suite@sha256:47eba58683641905a58f31c05605c953dc1d888466e8d434a6ceff47b76df03e"
        runtime_minutes: "600"
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

task AlleleCountsLocalize {
    input {
        String pairName
        String alleleCountsTxtPath = "~{pairName}.haplotypecaller.gatk.alleles.txt"
        IndexedReference referenceFa
        Bam normalFinalBam
        File knownHetVcf
        Bam tumorFinalBam
        Int memoryGb = 4

        Int diskSize
        
        Int timeLimit = 600
    }

    command {
    
        timeout ~{timeLimit}m \
        python3 \
        /parse_bam_generate_features.py \
        --tumor_bam ~{tumorFinalBam.bam} \
        --normal_bam ~{normalFinalBam.bam} \
        --vcf ~{knownHetVcf} \
        --output ~{alleleCountsTxtPath} \
        --reference ~{referenceFa.fasta}
    }

    output {
        File alleleCountsTxt = "~{alleleCountsTxtPath}"
    }

    runtime {
        mem: memoryGb + "G"
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-comp-s-fd4e/gcs_htslib_suite@sha256:47eba58683641905a58f31c05605c953dc1d888466e8d434a6ceff47b76df03e"
        runtime_minutes: "600"
    }
    
    parameter_meta {
        normalFinalBam: {
            description: "Input BAM filename",
            category: "required",
            localization_optional: false
        }
        tumorFinalBam: {
            description: "Input BAM filename",
            category: "required",
            localization_optional: false
        }
    }
}

task CalcBaf {
    input {
        String pairName
        String bafTxtPath = "~{pairName}.haplotypecaller.gatk.baf.txt"
        File alleleCountsTxt
        Int memoryGb = 4
        Int diskSize = (ceil( size(alleleCountsTxt, "GB") )  * 2 ) + 10
    }

    command {
        python \
        /calc_baf.py \
        ~{alleleCountsTxt} \
        ~{bafTxtPath}
    }

    output {
        File bafTxt = "~{bafTxtPath}"
    }

    runtime {
        mem: memoryGb + "G"
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/somatic_dna_tools@sha256:f281e73ddf515f3a5db6e766ef12fb2331dafa2b27c88e426764dcd8b4a7e19b"
        runtime_minutes: "90"
    }
}
