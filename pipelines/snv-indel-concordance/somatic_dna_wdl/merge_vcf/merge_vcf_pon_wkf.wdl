version 1.0

import "prep_merge_vcf_pon_wkf.wdl" as prepMergeVcfPon 
import "merge_callers_pon_wkf.wdl" as mergeCallersPon 
import "merge_vcf.wdl" as merge_vcf
import "../wdl_structs.wdl"


workflow MergeVcfPonWgs {
    # command
    #   Call variants in BAMs
    #   merge and filter raw VCFs
    #   annotate
    input {
        
        File filteredMantaSV 
        File mutect2
        String tumor
        
        IndexedReference referenceFa
        Array[String]+ listOfChroms
    }
    
    call prepMergeVcfPon.PrepMergeVcfPon as filteredMantaSVPrepMergeVcfPon {
        input:
            callerVcf=filteredMantaSV,
            tumor=tumor,
            tool='manta',
            referenceFa=referenceFa
            
    }
    
    call prepMergeVcfPon.PrepMergeVcfPon as mutect2PrepMergeVcfPon {
        input:
            callerVcf=mutect2,
            tumor=tumor,
            tool='mutect2',
            referenceFa=referenceFa
            
    }
    
    Array[IndexedVcf] allVcfCompressed = [filteredMantaSVPrepMergeVcfPon.preppedVcf, mutect2PrepMergeVcfPon.preppedVcf]
    

    call mergeCallersPon.MergeCallersPon {
        input:
            tumor = tumor,
            listOfChroms = listOfChroms,
            allVcfCompressed = allVcfCompressed
        
    }
    
    # in place of MergeChroms because no reorder step is needed
    call merge_vcf.Gatk4MergeSortVcf {
            input:
                tempVcfs = MergeCallersPon.finalChromVcf,
                sortedVcfPath = "~{tumor}.sorted.merged.v7.vcf",
                referenceFa = referenceFa,
                memoryGb = 8,
                diskSize = 10
    }
    
    output {
        File mergedVcf = Gatk4MergeSortVcf.sortedVcf.vcf
    }
}
