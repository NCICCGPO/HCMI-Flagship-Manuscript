version 1.0

import "calling.wdl" as calling
import "../wdl_structs.wdl"

workflow Manta {
    # command
    #   run Manta caller
    input {
        String tumor
        String normal
        String intHVmem = "unlimited"
        IndexedReference referenceFa
        IndexedTable callRegions
        Bam normalFinalBam
        Bam tumorFinalBam
        # reorder columns
        String pairName
        String mantaPath = "~{pairName}.manta.vcf"
        String filteredMantafPath = "~{pairName}.manta.filtered.vcf"
        # resources
        Int diskSize = ceil( size(tumorFinalBam.bam, "GB") + size(normalFinalBam.bam, "GB")) + 20
        Int memoryGb = 64
        Int threads = 8
        # remove definition after replacing the command step for gcp
        File mantaJsonLog
        
        Boolean highMem = false

    }
    
    Int lowCallMemoryGb = 4
    
    if (highMem) {
        Int highCallMemoryGb = 64
    }
    Int callMemoryGb = select_first([highCallMemoryGb, lowCallMemoryGb])

    call calling.MantaWgs {
        input:
            intHVmem = intHVmem,
            referenceFa = referenceFa,
            pairName = pairName,
            callRegions = callRegions,
            normalFinalBam = normalFinalBam,
            tumorFinalBam = tumorFinalBam,
            memoryGb = callMemoryGb,
            diskSize = diskSize,
            threads = threads
    }

    call calling.AddVcfCommand as mantaAddVcfCommand {
        input:
            inVcf = MantaWgs.somaticSV.vcf,
            jsonLog = mantaJsonLog,
            memoryGb = 4,
            diskSize = 10
    }

    call calling.ReorderVcfColumns as mantaReorderVcfColumns {
        input:
            tumor = tumor,
            normal = normal,
            rawVcf = mantaAddVcfCommand.outVcf,
            orderedVcfPath = mantaPath,
            memoryGb = 4,
            diskSize = 10
    }

    call calling.FilterNonpass {
        input:
            referenceFa = referenceFa,
            pairName = pairName,
            vcf = mantaReorderVcfColumns.orderedVcf,
            memoryGb = 8,
            threads = 4,
            diskSize = 10

    }

    call calling.ReorderVcfColumns as mantaFilteredReorderVcfColumns {
        input:
            tumor = tumor,
            normal = normal,
            rawVcf = FilterNonpass.outVcf,
            orderedVcfPath = filteredMantafPath,
            memoryGb = 4,
            diskSize = 10
    }

    output {
        IndexedVcf candidateSmallIndels = MantaWgs.candidateSmallIndels
        IndexedVcf diploidSV = MantaWgs.diploidSV
        IndexedVcf somaticSV = MantaWgs.somaticSV
        IndexedVcf candidateSV = MantaWgs.candidateSV
        File unfilteredMantaSV = mantaReorderVcfColumns.orderedVcf
        File filteredMantaSV = mantaFilteredReorderVcfColumns.orderedVcf
    }

}
