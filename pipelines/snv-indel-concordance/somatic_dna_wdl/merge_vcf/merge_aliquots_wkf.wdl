version 1.0

import "merge_vcf.wdl" as mergeVcf
import "../calling/calling.wdl" as calling
import "../calling/mutect2_multi_geno_wkf.wdl" as mutect2MultiGenotype
import "../tasks/utils.wdl" as utils
import "../wdl_structs.wdl"

# note that we will still need to match dictionaries or provide matching reference files

workflow MergeAliquots {
    input {
        Boolean local = false
        Array[IndexedVcf] centerVcfs
        Array[String] tumorBarcodeAliquots
        String normalBarcodeAliquot
        Array[Bam] tumorFinalBams
        Bam normalFinalBam
        String participantId
        File cosmicCensus
        
        Array[String]+ listOfChroms

        IndexedReference referenceFa
        Boolean reAnnotate = false
        
        File? serviceAccountKey
        String? gcpProject
    }

    if (length(centerVcfs) > 1 ) {
        scatter(vcfCompressed in centerVcfs) {
                File allVcfCompressedFile = vcfCompressed.vcf
            }
        Array[File]+ allVcfCompressedList = allVcfCompressedFile
        String infoTags = "aliquots_called_by:join,non_normal_kinds_called_by:join,MNV_ID:join"
        call mergeVcf.MergeCenterCallers {
                input:
                    pairName = participantId,
                    allVcfCompressed = centerVcfs,
                    allVcfCompressedList = allVcfCompressedList,
                    infoTags = infoTags
        }
    }
    File mergedVcf = select_first([MergeCenterCallers.mergedVcf, centerVcfs[0].vcf])
    # gather split VCF filenames to avoid glob that can fail on prem
    Array[String] suffixes = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29", "30", "31", "32", "33", "34", "35", "36", "37", "38", "39", "40", "41", "42", "43", "44", "45", "46", "47", "48", "49", "50"]
    Int maxSplits = 20
    String prefix = "~{participantId}.snv.indel.participant."
    String additionalSuffix = ".vcf"
    Array[String] mnvVcfPaths = ["~{prefix}.mnv~{additionalSuffix}"]
    scatter (index in range(length(suffixes))) {
        String suffix = suffixes[index]
        String splitChromVcfPaths = "~{prefix}.~{suffix}~{additionalSuffix}"
    }
    Array[String] splitVcfPaths = flatten([mnvVcfPaths, splitChromVcfPaths])
    call utils.SplitVcf {
        input:
            vcf=mergedVcf,
            prefix=prefix,
            diskSize = (ceil(size(mergedVcf, "GB")) * 3) + 10,
            maxRows = 1000,
            maxSplits = maxSplits,
            minSplits = 2,
            splitVcfPaths = splitVcfPaths
    }
    Array[File] splitVcfs = select_all(SplitVcf.splitVcfs)
    scatter(vcf in splitVcfs) {
        String splitId = sub(basename(vcf), ".vcf$", "")
        call mergeVcf.VcfToBed {
                input:
                    pairName = participantId,
                    chrom = splitId,
                    candidateChromVcf = vcf,
                    candidateChromBedPath = "~{participantId}.variant.bed"
            }
        scatter (i in range(length(tumorFinalBams))) {
            call calling.SliceBam as sliceBamTumor {
                    input:
                        chrom = splitId,
                        sampleId = tumorBarcodeAliquots[i],
                        chromBed = VcfToBed.candidateChromBed,
                        finalBam = tumorFinalBams[i],
                        threads = 4,
                        memoryGb = 4,
                        diskSize = 10,
                        serviceAccountKey = serviceAccountKey,
                        gcpProject = gcpProject
            }
            File slicedBamTumors = sliceBamTumor.lancetBamSlice.bam
        }
            
        call calling.SliceBam as sliceBamNormal {
            input:
                chrom = splitId,
                sampleId = normalBarcodeAliquot,
                chromBed = VcfToBed.candidateChromBed,
                finalBam = normalFinalBam,
                threads = 4,
                memoryGb = 4,
                diskSize = 10,
                serviceAccountKey = serviceAccountKey,
                gcpProject = gcpProject
        }
    
       Int diskSize = ceil( size(sliceBamNormal.lancetBamSlice.bam, "GB") + size(slicedBamTumors, "GB")) + 10
         call mergeVcf.CountAliquots {
             input:
                 splitId = splitId,
                 participantId = participantId,
                 mergedVcf = vcf,
                 normalFinalBam = sliceBamNormal.lancetBamSlice,
                 tumorFinalBams = sliceBamTumor.lancetBamSlice,
                 tumorFinalBamList = slicedBamTumors,
                 diskSize = diskSize
         }
    }
    
    Int mergeDiskSize = (ceil(size(CountAliquots.countedVcf, "GB")) * 2) + 10
    call mergeVcf.MergeBcftools {
        input:
            vcfs = CountAliquots.countedVcf,
            sortedVcfPath = "~{participantId}.snv.indel.counted.vcf",
            memoryGb = 16,
            diskSize = mergeDiskSize
    }
    
    call mutect2MultiGenotype.Mutect2MultiGenotype {
            input:
                local = local,
                tumors = tumorBarcodeAliquots,
                normal = normalBarcodeAliquot,
                allelesVcf = mergedVcf,
                participantId = participantId,
                referenceFa = referenceFa,
                normalFinalBam = normalFinalBam,
                tumorFinalBams = tumorFinalBams,
                
    }
    
    Int diskSizeWithMutect2 = (ceil( size(MergeBcftools.sortedVcf, "GB") + size(Mutect2MultiGenotype.mutect2, "GB") )  * 4 ) + 50
    call mergeVcf.AnnotateWithMutect2 {
        input:
            mergedVcf = MergeBcftools.sortedVcf,
            mutect2 = Mutect2MultiGenotype.mutect2,
            participantId = participantId,
            cosmicCensus = cosmicCensus,
            diskSize = diskSizeWithMutect2
    }
    
    output {
        File unionVcf = AnnotateWithMutect2.unionVcf
        File intersectionVcf = AnnotateWithMutect2.intersectionVcf
        File mutect2Genotyped = Mutect2MultiGenotype.mutect2
        File mutect2UnfilteredGenotyped = Mutect2MultiGenotype.mutect2Unfiltered
        }
}
