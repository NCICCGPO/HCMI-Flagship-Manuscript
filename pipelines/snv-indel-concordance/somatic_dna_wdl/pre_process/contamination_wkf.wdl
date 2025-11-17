version 1.0

import "../wdl_structs.wdl"
import "qc.wdl"

workflow CalculateContamination {
    # command
    input {
        Bam finalNormalBam
        Bam finalTumorBam
        IndexedReference referenceFa
        String tumor
        String normal
        String pairName

        IndexedVcf gnomadBiallelic
        Int threads = 1
        Int memoryGb = 16
    }

    Int additionalDiskSize = 100
    Int tumorSize = ceil(size(finalTumorBam.bam, "GB") + size(finalTumorBam.bamIndex, "GB"))
    Int normalSize = ceil(size(finalNormalBam.bam, "GB") + size(finalNormalBam.bamIndex, "GB"))

    call qc.Pileup as normalPileup {
        input:
            sampleId = normal,
            finalBam = finalNormalBam,
            gnomadBiallelic = gnomadBiallelic,
            memoryGb = memoryGb,
            diskSize = normalSize + additionalDiskSize
    }

    call qc.Pileup as tumorPileup {
        input:
            sampleId = tumor,
            finalBam = finalTumorBam,
            gnomadBiallelic = gnomadBiallelic,
            memoryGb = memoryGb,
            diskSize = tumorSize + additionalDiskSize
    }

    call qc.CalculateContaminationPaired {
        input:
            pairName = pairName,
            pileupsNormalTable = normalPileup.pileupsTable,
            pileupsTumorTable = tumorPileup.pileupsTable,
            memoryGb = memoryGb,
            threads = threads
    }
}
