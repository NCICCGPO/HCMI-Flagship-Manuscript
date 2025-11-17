version 1.0

import "calling.wdl" as calling
import "../wdl_structs.wdl"

workflow MantaPon {
    # command
    #   run Manta caller
    input {
        String tumor
        String intHVmem = "unlimited"
        IndexedReference referenceFa
        IndexedTable callRegions
        Bam tumorFinalBam
        # reorder columns
        String mantaPath = "~{tumor}.manta.vcf"
        String filteredMantafPath = "~{tumor}.manta.filtered.vcf"
        # resources
        Int diskSize = ceil( size(tumorFinalBam.bam, "GB")) + 20
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

    call calling.MantaWgsPon {
        input:
            intHVmem = intHVmem,
            referenceFa = referenceFa,
            sampleId = tumor,
            callRegions = callRegions,
            tumorFinalBam = tumorFinalBam,
            memoryGb = callMemoryGb,
            diskSize = diskSize,
            threads = threads
    }
    
    call calling.FilterNonpassPon {
        input:
            referenceFa = referenceFa,
            pairName = tumor,
            vcf = MantaWgsPon.tumorSV,
            memoryGb = 8,
            outVcfPath = "~{tumor}.manta.filtered.vcf",
            threads = 4,
            diskSize = 10

    }

    output {
        File filteredMantaSV = FilterNonpassPon.outVcf
    }

}
