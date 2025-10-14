version 1.0

import "merge_vcf.wdl" as merge_vcf
import "../calling/calling.wdl" as calling
import "../wdl_structs.wdl"

workflow MergeChroms {
    input {
        String tumorId
        String normalId
        String pairName
        String sortedVcfPath = "~{pairName}.sorted.merged.v7.vcf"
        String orderedVcfPath = "~{pairName}.merged.v7.vcf"
        Array[File]+ finalChromVcf
        IndexedReference referenceFa
    }

    call merge_vcf.Gatk4MergeSortVcf {
            input:
                tempVcfs = finalChromVcf,
                sortedVcfPath = sortedVcfPath,
                referenceFa = referenceFa,
                memoryGb = 8,
                diskSize = 10
        }

    call calling.ReorderVcfColumns {
        input:
            tumor = tumorId,
            normal = normalId,
            rawVcf = Gatk4MergeSortVcf.sortedVcf.vcf,
            orderedVcfPath = orderedVcfPath,
            memoryGb = 4,
            diskSize = 10
    }

    output {
        File unannotatedVcf = ReorderVcfColumns.orderedVcf
    }

}
