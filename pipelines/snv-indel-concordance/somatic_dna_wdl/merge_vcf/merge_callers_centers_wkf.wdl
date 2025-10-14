version 1.0

import "merge_vcf.wdl" as mergeVcf
import "../calling/calling.wdl" as calling
import "../tasks/utils.wdl" as utils
import "../wdl_structs.wdl"

# note that we will still need to match dictionaries or provide matching reference files

workflow MergeCenterCallers {
    input {
        Boolean external = false
        String tumorId
        String normalId
        String pairName
        Array[String]+ listOfChroms
        Array[IndexedVcf]+ allVcfCompressed

        File intervalListBed

        IndexedReference referenceFa
        Bam normalFinalBam
        Bam tumorFinalBam
        
        String infoTags
        File? serviceAccountKey
        String? gcpProject
    }

    scatter(vcfCompressed in allVcfCompressed) {
            File allVcfCompressedFile = vcfCompressed.vcf
        }
    Array[File]+ allVcfCompressedList = allVcfCompressedFile
    
    call mergeVcf.MergeCenterCallers {
            input:
                pairName = pairName,
                allVcfCompressed = allVcfCompressed,
                allVcfCompressedList = allVcfCompressedList,
                infoTags = infoTags
    }
    
    # gather split VCF filenames to avoid glob that can fail on prem
    Array[String] suffixes = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29", "30", "31", "32", "33", "34", "35", "36", "37", "38", "39", "40", "41", "42", "43", "44", "45", "46", "47", "48", "49", "50"]
    Int maxSplits = 20
    String prefix = "~{pairName}.snv.indel.v7."
    String additionalSuffix = ".vcf"
    Array[String] mnvVcfPaths = ["~{prefix}.mnv~{additionalSuffix}"]
    scatter (index in range(length(suffixes))) {
        String suffix = suffixes[index]
        String splitChromVcfPaths = "~{prefix}.~{suffix}~{additionalSuffix}"
    }
    Array[String] splitVcfPaths = flatten([mnvVcfPaths, splitChromVcfPaths])
    call utils.SplitVcf {
        input:
            vcf=MergeCenterCallers.mergedVcf,
            prefix=prefix,
            diskSize = (ceil(size(MergeCenterCallers.mergedVcf, "GB")) * 3) + 10,
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
                    pairName = pairName,
                    chrom = splitId,
                    candidateChromVcf = vcf,
                    candidateChromBedPath = "~{pairName}.variant.bed"
            }
        call calling.SliceBam as sliceBamTumor {
                input:
                    chrom = splitId,
                    sampleId = tumorId,
                    chromBed = VcfToBed.candidateChromBed,
                    finalBam = tumorFinalBam,
                    threads = 4,
                    memoryGb = 4,
                    diskSize = 10,
                    serviceAccountKey = serviceAccountKey,
                    gcpProject = gcpProject
            }
        
        call calling.SliceBam as sliceBamNormal {
            input:
                chrom = splitId,
                sampleId = normalId,
                chromBed = VcfToBed.candidateChromBed,
                finalBam = normalFinalBam,
                threads = 4,
                memoryGb = 4,
                diskSize = 10,
                serviceAccountKey = serviceAccountKey,
                gcpProject = gcpProject
        }
    
       Int diskSize = ceil( size(sliceBamNormal.lancetBamSlice.bam, "GB") + size(sliceBamTumor.lancetBamSlice.bam, "GB")) + 10
  
        call mergeVcf.CountCenterCallers {
            input:
                chrom = splitId,
                pairName = pairName,
                mergedVcf = vcf,
                tumorId = tumorId,
                normalId = normalId,
                normalFinalBam = sliceBamNormal.lancetBamSlice,
                tumorFinalBam = sliceBamTumor.lancetBamSlice,
                diskSize = diskSize
        }
    }

    output {
            Array[File] mergedChromVcf = CountCenterCallers.mergedChromVcf
            Array[File] finalChromVcf = CountCenterCallers.finalChromVcf
        }
}
