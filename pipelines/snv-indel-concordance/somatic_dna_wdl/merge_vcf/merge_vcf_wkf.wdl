version 1.0

import "prep_merge_vcf_wkf.wdl" as prepMergeVcf
import "merge_callers_wkf.wdl" as mergeCallers
import "merge_chroms_wkf.wdl" as mergeChroms
import "../tasks/utils.wdl" as utils
import "../wdl_structs.wdl"

workflow MergeVcf {
    # command
    #   Call variants in BAMs
    #   merge and filter raw VCFs
    #   annotate
    input {
        Boolean external = false
        String library
        PreMergedPairVcfInfo preMergedPairVcfInfo
        IndexedReference referenceFa
        Array[String]+ listOfChroms
        
        # merge callers
        File intervalListBed
        File ponFile
        
        # needed for work on cloud without localizing
        File? serviceAccountKey
        String? gcpProject
    }
    
    if (library == 'WGS') {
        call prepMergeVcf.PrepMergeVcf as filteredMantaSVPrepMergeVcf {
            input:
                callerVcf=preMergedPairVcfInfo.filteredMantaSV,
                tumor=preMergedPairVcfInfo.tumor,
                normal=preMergedPairVcfInfo.normal,
                tool='manta',
                pairName=preMergedPairVcfInfo.pairId,
                referenceFa=referenceFa
                
        }
    }
    
    if (library == 'Exome') {
        call utils.CreateBlankFile as createVcf {
                input:
                fileId = "~{preMergedPairVcfInfo.pairId}_nonMantaSvVcf_"
        }
        
        call utils.CreateBlankFile as createVcfIndex {
            input:
                fileId = "~{preMergedPairVcfInfo.pairId}_nonMantaSvVcfIndex_"
        }
        
        IndexedVcf nonfilteredMantaSV = object {
                vcf : createVcf.blankFile,
                index : createVcfIndex.blankFile
        }
    }
    
    IndexedVcf preppedFilteredMantaSV = select_first([filteredMantaSVPrepMergeVcf.preppedVcf, nonfilteredMantaSV])
    
    call prepMergeVcf.PrepMergeVcf as strelka2SnvPrepMergeVcf {
        input:
            callerVcf=preMergedPairVcfInfo.strelka2Snv,
            tumor=preMergedPairVcfInfo.tumor,
            normal=preMergedPairVcfInfo.normal,
            tool='strelka2',
            pairName=preMergedPairVcfInfo.pairId,
            referenceFa=referenceFa
            
    }
    
    call prepMergeVcf.PrepMergeVcf as strelka2IndelPrepMergeVcf {
        input:
            callerVcf=preMergedPairVcfInfo.strelka2Indel,
            tumor=preMergedPairVcfInfo.tumor,
            normal=preMergedPairVcfInfo.normal,
            tool='strelka2',
            pairName=preMergedPairVcfInfo.pairId,
            referenceFa=referenceFa
            
    }
    
    call prepMergeVcf.PrepMergeVcf as mutect2PrepMergeVcf {
        input:
            callerVcf=preMergedPairVcfInfo.mutect2,
            tumor=preMergedPairVcfInfo.tumor,
            normal=preMergedPairVcfInfo.normal,
            tool='mutect2',
            pairName=preMergedPairVcfInfo.pairId,
            referenceFa=referenceFa
            
    }
    
    call prepMergeVcf.PrepMergeVcf as lancetPrepMergeVcf {
        input:
            callerVcf=preMergedPairVcfInfo.lancet,
            tumor=preMergedPairVcfInfo.tumor,
            normal=preMergedPairVcfInfo.normal,
            tool='lancet',
            pairName=preMergedPairVcfInfo.pairId,
            referenceFa=referenceFa
            
    }
    
    if (library == 'WGS') {
        Array[IndexedVcf] allVcfCompressedWgs = [preppedFilteredMantaSV, 
                    strelka2SnvPrepMergeVcf.preppedVcf,
                    strelka2IndelPrepMergeVcf.preppedVcf,
                    mutect2PrepMergeVcf.preppedVcf,
                    lancetPrepMergeVcf.preppedVcf]
    }
    
    if (library == 'Exome') {
        Array[IndexedVcf] allVcfCompressedExome = [strelka2SnvPrepMergeVcf.preppedVcf,
                    strelka2IndelPrepMergeVcf.preppedVcf,
                    mutect2PrepMergeVcf.preppedVcf,
                    lancetPrepMergeVcf.preppedVcf]
    }
    
    Array[IndexedVcf] allVcfCompressed = select_first([allVcfCompressedWgs, allVcfCompressedExome])
    
    call mergeCallers.MergeCallers {
        input:
            external=external,
            tumorId=preMergedPairVcfInfo.tumor,
            normalId=preMergedPairVcfInfo.normal,
            pairName=preMergedPairVcfInfo.pairId,
            listOfChroms=listOfChroms,
            intervalListBed=intervalListBed,
            referenceFa=referenceFa,
            normalFinalBam=preMergedPairVcfInfo.normalFinalBam,
            tumorFinalBam=preMergedPairVcfInfo.tumorFinalBam,
            ponFile=ponFile,
            allVcfCompressed=allVcfCompressed,
            serviceAccountKey = serviceAccountKey,
            gcpProject = gcpProject
        
    }
    
    call mergeChroms.MergeChroms {
        input:
            tumorId=preMergedPairVcfInfo.tumor,
            normalId=preMergedPairVcfInfo.normal,
            pairName=preMergedPairVcfInfo.pairId,
            referenceFa=referenceFa,
            finalChromVcf=MergeCallers.finalChromVcf,
            
    }
    
    output {
        File mergedVcf = MergeChroms.unannotatedVcf
    }
}
