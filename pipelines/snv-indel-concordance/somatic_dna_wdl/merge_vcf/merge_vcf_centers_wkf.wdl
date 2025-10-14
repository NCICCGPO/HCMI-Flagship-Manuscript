version 1.0

import "prep_merge_vcf_centers_wkf.wdl" as prepMergeVcf
import "merge_callers_centers_wkf.wdl" as mergeCallers
import "merge_chroms_wkf.wdl" as mergeChroms
import "merge_vcf.wdl" as mergeVcf
import "../calling/calling.wdl" as calling
import "../wdl_structs.wdl"

workflow MergeVcfCenters {
    # command
    #   Call variants in BAMs
    #   merge and filter raw VCFs
    #   annotate
    input {
        Boolean external = false
        File? broadMaf
        File? nygcVcf
        File? washuVcf
        String tumorId
        String normalId
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
    
    if (defined(broadMaf)) {
        # adjust to avoid failure for pipeline failures
        File runBroadMaf = select_first([broadMaf, referenceFa.fasta])
        call mergeVcf.PrepBroadMaf {
            input:
                broadMaf = runBroadMaf,
                referenceFa = referenceFa.fasta,
                tumor = tumorId,
                normal = normalId
        }
        call mergeVcf.RenameBroadVcf {
            input:
                cleanVcf = PrepBroadMaf.broadCleanVcf,
                tumor = tumorId,
                normal = normalId
        }
        call mergeVcf.Gatk4MergeSortVcf as broadGatk4MergeSortVcf {
            input:
                tempVcfs = [RenameBroadVcf.broadVcf],
                sortedVcfPath = sub(basename(RenameBroadVcf.broadVcf), ".vcf$", ".contig.vcf"),
                referenceFa = referenceFa,
                gzipped = false,
                threads = 4,
                memoryGb = 8,
                diskSize = 10
    
        }
        
        call calling.ReorderVcfColumns {
            input:
                tumor = tumorId,
                normal = normalId,
                rawVcf = broadGatk4MergeSortVcf.sortedVcf.vcf,
                orderedVcfPath = sub(basename(broadGatk4MergeSortVcf.sortedVcf.vcf), ".vcf$", ".contig.reordered.vcf"),
                memoryGb = 4,
                diskSize = 10
        }    
        call prepMergeVcf.PrepMergeVcf as broadPrepMergeVcf {
            input:
                callerVcf=ReorderVcfColumns.orderedVcf,
                tumorId=tumorId,
                normalId=normalId,
                tool='broad',
                pairName=pairId,
                referenceFa=referenceFa,
                library=library
    
        }
    }
    
    # adjust to avoid failure for pipeline failures
    File runNygcVcf = select_first([nygcVcf, referenceFa.fasta])
    if (defined(nygcVcf)) {
        call mergeVcf.PrepNygcVcf {
            input:
                nygcVcf = runNygcVcf,
                normals = [normalId]
        }
        call prepMergeVcf.PrepMergeVcf as nygcPrepMergeVcf {
            input:
                callerVcf=PrepNygcVcf.nygcPreppedVcf,
                tumorId=tumorId,
                normalId=normalId,
                tool='nygc',
                pairName=pairId,
                referenceFa=referenceFa,
                library=library
                
        }
    }
    # adjust to avoid failure for pipeline failures
    File runWashuVcf = select_first([washuVcf, referenceFa.fasta])
    if (defined(washuVcf)) {
        call mergeVcf.PrepWashuVcf {
            input:
                washuVcf = runWashuVcf,
                tumor = tumorId,
                normal = normalId
        }
        call prepMergeVcf.PrepMergeVcf as washuPrepMergeVcf {
            input:
                callerVcf=PrepWashuVcf.washuPreppedVcf,
                tumorId=tumorId,
                normalId=normalId,
                tool='washu',
                pairName=pairId,
                referenceFa=referenceFa,
                library=library
                
        }
    }
    
    if (library == "WGS"){
        String wgsInfoTags = "centers_called_by:join,num_centers:sum,called_by:join,MNV_ID:join"
    }
    if (library == "WgsExome"){
        String exomeInfoTags = "centers_called_by:join,num_centers:sum,Exome_called_by:join,WGS_called_by:sum,called_by:join,num_callers:sum,MNV_ID:join"
    }
    
    if (library == "Exome"){
        String exomeOnlyInfoTags = "centers_called_by:join,num_centers:sum,called_by:join,MNV_ID:join"
    }
    String infoTags = select_first([wgsInfoTags, exomeInfoTags, exomeOnlyInfoTags])
    
    Array[IndexedVcf] allVcfCompressed = select_all([washuPrepMergeVcf.preppedVcf, nygcPrepMergeVcf.preppedVcf, broadPrepMergeVcf.preppedVcf])

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
            allVcfCompressed=allVcfCompressed,
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
