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
        String library
        IndexedReference referenceFa
    }
    
    if (library == "WGS") {
        call merge_vcf.PrepCenterCalls {
            input:
                callerVcf = callerVcf,
                tool=tool,
                normalId=normalId,
                tumorId=tumorId,
                pairName=pairName,
                referenceFa=referenceFa
        }
    }
    if (library == "WgsExome") {
        call merge_vcf.PrepCenterCallsExomeWgs {
            input:
                callerVcf = callerVcf,
                tool=tool,
                normalId=normalId,
                tumorId=tumorId,
                pairName=pairName,
                referenceFa=referenceFa
        }
    }
    
    File mnvVcf = select_first([PrepCenterCalls.mnvVcf, PrepCenterCallsExomeWgs.mnvVcf])
    
    call merge_vcf.Gatk4MergeSortVcf {
        input:
            tempVcfs = [mnvVcf],
            sortedVcfPath = sub(basename(mnvVcf), "$", ".gz"),
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
