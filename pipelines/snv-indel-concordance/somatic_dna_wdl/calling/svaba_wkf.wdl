version 1.0

import "calling.wdl" as calling
import "../wdl_structs.wdl"

workflow Svaba {
    # command
    #   run Svaba caller
    input {
        String tumor
        String normal
        String pairName
        BwaReference bwaReference
        IndexedTable callRegions
        File dbsnpIndels
        Bam normalFinalBam
        Bam tumorFinalBam
        # resources
        Int threads = 8
        Int diskSize = ceil( size(tumorFinalBam.bam, "GB") + size(normalFinalBam.bam, "GB")) + 40
        Int memoryGb = 32
        File svabaJsonLog
        
        Boolean highMem = false
    }
    
    Int lowCallMemoryGb = 16
    
    if (highMem) {
        Int highCallMemoryGb = 32
    }
    Int callMemoryGb = select_first([highCallMemoryGb, lowCallMemoryGb])

    call calling.SvabaWgs {
        input:
            pairName=pairName,
            bwaReference = bwaReference,
            callRegions = callRegions,
            normalFinalBam = normalFinalBam,
            tumorFinalBam = tumorFinalBam,
            dbsnpIndels = dbsnpIndels,
            memoryGb = callMemoryGb,
            threads = threads,
            diskSize = diskSize

    }

    # sv
    call calling.AddVcfCommand as svAddVcfCommand {
        input:
            inVcf = SvabaWgs.svabaGz,
            jsonLog = svabaJsonLog,
            memoryGb = 4,
            diskSize = 10
    }

    call calling.ReorderVcfColumns as svReorderVcfColumns {
        input:
            tumor = tumor,
            normal = normal,
            rawVcf = svAddVcfCommand.outVcf,
            orderedVcfPath = "~{pairName}.sv.svaba.vcf",
            memoryGb = 4,
            diskSize = 10
    }

    # indel
    call calling.AddVcfCommand as indelAddVcfCommand {
        input:
            inVcf = SvabaWgs.svabaIndelGz,
            jsonLog = svabaJsonLog,
            memoryGb = 4,
            diskSize = 10
    }

    call calling.ReorderVcfColumns as indelReorderVcfColumns {
        input:
            tumor = tumor,
            normal = normal,
            rawVcf = indelAddVcfCommand.outVcf,
            orderedVcfPath = "~{pairName}.indel.svaba.vcf",
            memoryGb = 4,
            diskSize = 10
    }

    output {
        File svabaRawGermlineIndel = SvabaWgs.svabaRawGermlineIndel
        File svabaRawGermlineSv = SvabaWgs.svabaRawGermlineSv
        File svabaSv = svReorderVcfColumns.orderedVcf
        File svabaIndel = indelReorderVcfColumns.orderedVcf
    }
}
