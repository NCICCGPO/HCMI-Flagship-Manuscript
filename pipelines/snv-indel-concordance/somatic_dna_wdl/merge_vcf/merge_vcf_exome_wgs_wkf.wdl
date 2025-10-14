version 1.0

import "prep_merge_vcf_exome_wgs_wkf.wdl" as prepMergeVcf
import "merge_callers_exome_wgs_wkf.wdl" as mergeCallers
import "merge_chroms_wkf.wdl" as mergeChroms
import "../wdl_structs.wdl"

workflow MergeVcfExomeWgs {
    # command
    #   Call variants in BAMs
    #   merge and filter raw VCFs
    #   annotate
    input {
        Boolean external = false
        File wgsVcf
        File exomeVcf
        String tumorId
        String normalId
        String pairId
        String exomeTumorId
        String exomeNormalId
        String exomePairId
        
        IndexedReference referenceFa
        Array[String]+ listOfChroms
    }
    
    call prepMergeVcf.PrepMergeVcf as wgsPrepMergeVcf {
        input:
            callerVcf=wgsVcf,
            tumorId=tumorId,
            normalId=normalId,
            tool='WGS',
            pairName=pairId,
            referenceFa=referenceFa
            
    }
    
    call prepMergeVcf.PrepMergeVcf as exomePrepMergeVcf {
        input:
            callerVcf=exomeVcf,
            tumorId=tumorId,
            normalId=normalId,
            tool='Exome',
            pairName=pairId,
            referenceFa=referenceFa
            
    }
    
    call mergeCallers.MergeCallers {
        input:
            external=external,
            tumorId=tumorId,
            normalId=normalId,
            pairName=pairId,
            listOfChroms=listOfChroms,
            referenceFa=referenceFa,
            allVcfCompressed=[wgsPrepMergeVcf.preppedVcf, exomePrepMergeVcf.preppedVcf]
        
    }
    
    call mergeChroms.MergeChroms {
        input:
            tumorId=tumorId,
            normalId=normalId,
            pairName=pairId,
            referenceFa=referenceFa,
            finalChromVcf=MergeCallers.finalChromVcf,
            
    }
    
    output {
        File mergedVcf = MergeChroms.unannotatedVcf
    }
}
