version 1.0


import "pon_merge.wdl" as ponMerge
import "../wdl_structs.wdl"

workflow PrepSinglePon {
    input {
        File vcfAnnotatedVep
        String sampleId
        Int diskSize
    }
    
    call ponMerge.MakeSinglePon {
        input:
            diskSize = diskSize,
            sampleId = sampleId,
            vcfAnnotatedVep = vcfAnnotatedVep
    }
    
    output {
        File singlePon = MakeSinglePon.singlePon
    }
}