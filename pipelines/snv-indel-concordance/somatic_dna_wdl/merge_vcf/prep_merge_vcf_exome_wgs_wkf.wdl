version 1.0

import "merge_vcf.wdl" as merge_vcf
import "../wdl_structs.wdl"

workflow PrepMergeVcf {
    input {
        String tumorId
        String normalId
        String tool
        File callerVcf
        String pairName
        IndexedReference referenceFa
    }
    
    call merge_vcf.PrepExomeWgsMetadata {
        input:
            callerVcf = callerVcf,
            tumorId = tumorId,
            normalId = normalId,
            referenceFa = referenceFa,
            tool=tool,
            pairName=pairName
    }
    
    call merge_vcf.Gatk4MergeSortVcf {
        input:
            tempVcfs = [PrepExomeWgsMetadata.mnvVcf],
            sortedVcfPath = sub(basename(PrepExomeWgsMetadata.mnvVcf), "$", ".gz"),
            referenceFa = referenceFa,
            gzipped = true,
            threads = 4,
            memoryGb = 8,
            diskSize = 20
            
    }
    
    output {
        IndexedVcf preppedVcf = Gatk4MergeSortVcf.sortedVcf
    }
}
