version 1.0

import "prep_merge_vcf_centers_wkf.wdl" as prepMergeVcf
import "merge_callers_centers_wkf.wdl" as mergeCallers
import "merge_chroms_wkf.wdl" as mergeChroms
import "merge_vcf.wdl" as mergeVcf
import "../wdl_structs.wdl"

workflow MergeVcfCenters {
    # command
    #   Call variants in BAMs
    #   merge and filter raw VCFs
    #   annotate
    input {
        Boolean external = false
        File broadVcf
        File broadExomeVcf
        File nygcVcf
        File washuVcf
        File washuExomeVcf
        String tumorId
        String normalId
        String tumorIdExome
        String normalIdExome
        String pairId
        
        Bam normalFinalBam
        Bam tumorFinalBam
        
        IndexedReference referenceFa
        Array[String]+ listOfChroms
        
        # merge callers
        File intervalListBed
        String library = "WGS"
        File? serviceAccountKey
        String? gcpProject
    }
    
    call mergeVcf.Gatk4MergeSortVcf as broadGatk4MergeSortVcf {
        input:
            tempVcfs = [broadVcf],
            sortedVcfPath = sub(basename(broadVcf), ".vcf$", ".contig.vcf"),
            referenceFa = referenceFa,
            gzipped = false,
            threads = 4,
            memoryGb = 8,
            diskSize = 10

    }
    call mergeVcf.Gatk4MergeSortVcf as broadExomeGatk4MergeSortVcf {
        input:
            tempVcfs = [broadExomeVcf],
            sortedVcfPath = sub(basename(broadExomeVcf), ".vcf$", ".contig.vcf"),
            referenceFa = referenceFa,
            gzipped = false,
            threads = 4,
            memoryGb = 8,
            diskSize = 10

    }
    
    call prepMergeVcf.PrepMergeVcf as broadPrepMergeVcf {
        input:
            callerVcf=broadGatk4MergeSortVcf.sortedVcf.vcf,
            tumorId=tumorId,
            normalId=normalId,
            tool='broadWgs',
            pairName=pairId,
            library=library,
            referenceFa=referenceFa

    }
    call prepMergeVcf.PrepMergeVcf as broadExomePrepMergeVcf {
        input:
            callerVcf=broadExomeGatk4MergeSortVcf.sortedVcf.vcf,
            tumorId=tumorId,
            normalId=normalId,
            tool='broadExome',
            pairName=pairId,
            library=library,
            referenceFa=referenceFa

    }

    call prepMergeVcf.PrepMergeVcf as nygcPrepMergeVcf {
        input:
            callerVcf=nygcVcf,
            tumorId=tumorId,
            normalId=normalId,
            tool='nygc',
            pairName=pairId,
            library=library,
            referenceFa=referenceFa
            
    }
    
    call mergeVcf.Gatk4MergeSortVcf as washuGatk4MergeSortVcf {
        input:
            tempVcfs = [washuVcf],
            sortedVcfPath = sub(basename(washuVcf), ".vcf$", ".contig.vcf"),
            referenceFa = referenceFa,
            gzipped = false,
            threads = 4,
            memoryGb = 8,
            diskSize = 10

    }
    call mergeVcf.Gatk4MergeSortVcf as washuExomeGatk4MergeSortVcf {
        input:
            tempVcfs = [washuExomeVcf],
            sortedVcfPath = sub(basename(washuExomeVcf), ".vcf$", ".contig.vcf"),
            referenceFa = referenceFa,
            gzipped = false,
            threads = 4,
            memoryGb = 8,
            diskSize = 10

    }
    call prepMergeVcf.PrepMergeVcf as washuPrepMergeVcf {
        input:
            callerVcf=washuGatk4MergeSortVcf.sortedVcf.vcf,
            tumorId=tumorId,
            normalId=normalId,
            tool='washuWgs',
            pairName=pairId,
            library=library,
            referenceFa=referenceFa
            
    }
    call prepMergeVcf.PrepMergeVcf as washuExomePrepMergeVcf {
        input:
            callerVcf=washuExomeGatk4MergeSortVcf.sortedVcf.vcf,
            tumorId=tumorId,
            normalId=normalId,
            tool='washuExome',
            pairName=pairId,
            library=library,
            referenceFa=referenceFa
            
    }
    
    if (library == "WGS"){
        String wgsInfoTags = "centers_called_by:join,num_centers:sum,called_by:join,MNV_ID:join"
    }
    if (library == "WgsExome"){
        String exomeInfoTags = "centers_called_by:join,num_centers:sum,Exome_called_by:join,WGS_called_by:join,called_by:join,MNV_ID:join"
    }
    String infoTags = select_first([wgsInfoTags, exomeInfoTags])
    call mergeCallers.MergeCenterCallers {
        input:
            external=external,
            tumorId=tumorId,
            normalId=normalId,
            pairName=pairId,
            infoTags=infoTags,
            listOfChroms=listOfChroms,
            intervalListBed=intervalListBed,
            referenceFa=referenceFa,
            normalFinalBam=normalFinalBam,
            tumorFinalBam=tumorFinalBam,
            allVcfCompressed=[broadPrepMergeVcf.preppedVcf, broadExomePrepMergeVcf.preppedVcf , nygcPrepMergeVcf.preppedVcf, washuPrepMergeVcf.preppedVcf, washuExomePrepMergeVcf.preppedVcf],
            serviceAccountKey = serviceAccountKey,
            gcpProject = gcpProject
    }
    
    call mergeChroms.MergeChroms {
        input:
            tumorId=tumorId,
            normalId=normalId,
            pairName=pairId,
            referenceFa=referenceFa,
            finalChromVcf=MergeCenterCallers.finalChromVcf,
    }
    
    output {
        File mergedVcf = MergeChroms.unannotatedVcf
    }
}
