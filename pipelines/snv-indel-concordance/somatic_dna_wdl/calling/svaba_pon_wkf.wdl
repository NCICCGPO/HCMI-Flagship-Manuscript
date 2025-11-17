version 1.0

import "calling.wdl" as calling
import "../wdl_structs.wdl"

workflow SvabaPon {
    # command
    #   run Svaba caller
    input {
        String tumor
        File refCache
        BwaReference svabaIndexedReference
        IndexedTable callRegions
        IndexedReference referenceFa
        File dbsnpIndels
        Bam tumorFinalBam
        # resources
        Int threads = 4
        Int diskSize = ceil( size(tumorFinalBam.bam, "GB")) + 20
        File svabaJsonLog
        
        Boolean highMem = true
    }
    
    Int lowCallMemoryGb = 16
    
    if (highMem) {
        Int highCallMemoryGb = 32
    }
    Int callMemoryGb = select_first([highCallMemoryGb, lowCallMemoryGb])


    call calling.SvabaWgsPon {
        input:
            refCache = refCache,
            sampleId=tumor,
            svabaIndexedReference = svabaIndexedReference,
            callRegions = callRegions,
            tumorFinalBam = tumorFinalBam,
            dbsnpIndels = dbsnpIndels,
            memoryGb = callMemoryGb,
            threads = threads,
            diskSize = diskSize
    }

    # indel
    call calling.ReheaderVcf as indelReheaderVcf {
        input:
            inVcf = SvabaWgsPon.svabaIndelGz,
            sampleIds = [tumor],
            referenceFa = referenceFa,
            outVcfPath = "~{tumor}.indel.svaba.vcf",
            memoryGb = 4,
            diskSize = 10
    }

    output {
        File svabaSvGz = SvabaWgsPon.svabaGz
        File svabaIndel = indelReheaderVcf.reheaderedVcf
    }
}
