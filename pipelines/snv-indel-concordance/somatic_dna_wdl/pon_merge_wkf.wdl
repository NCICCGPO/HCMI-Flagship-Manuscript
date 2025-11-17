version 1.0


import "merge_bed/pon_merge_wkf.wdl" as makeSinglePon
import "merge_bed/pon_merge.wdl" as ponMergeTasks
import "wdl_structs.wdl"

workflow PonMerge {
    input {
        Array[File] vcfAnnotatedVeps
        Array[String] sampleIds
    }
    
    scatter (vcfAnnotatedVep in zip(vcfAnnotatedVeps, sampleIds)) { 
        Int additionalDiskSize = 5
        Int diskSize = ceil((size(vcfAnnotatedVep.left, "GB"))) + additionalDiskSize
               
        call makeSinglePon.PrepSinglePon {
            input:
                sampleId = vcfAnnotatedVep.right,
                vcfAnnotatedVep = vcfAnnotatedVep.left,
                diskSize = diskSize
        }
    }
    
    call ponMergeTasks.MergePons {
        input:
            diskSize = 300,
            singlePons = PrepSinglePon.singlePon
    }
    
    output {
        File pon = MergePons.pon
    }
}