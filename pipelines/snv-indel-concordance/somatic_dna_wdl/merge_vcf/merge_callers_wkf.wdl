version 1.0

import "merge_vcf.wdl" as mergeVcf
import "../calling/calling.wdl" as calling
import "../tasks/utils.wdl" as utils
import "../wdl_structs.wdl"

# note that we will still need to match dictionaries or provide matching reference files

workflow MergeCallers {
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

        File ponFile

        Int lancetSplits = 30
        # needed for work on cloud without localizing
        File? serviceAccountKey
        String? gcpProject
    }
    
    # gather split BED filenames to avoid glob that can fail on prem
    Array[String] suffixes = ["00", "01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29", "30", "31", "32", "33", "34", "35", "36", "37", "38", "39", "40", "41", "42", "43", "44", "45", "46", "47", "48", "49", "50"]
    String candidateBedPrefixPath = "~{pairName}.candidate.split."
    String additionalSuffix = ".bed"
    scatter (index in range(length(suffixes))) {
        String suffix = suffixes[index]
        String candidateBedPaths = "~{candidateBedPrefixPath}~{suffix}~{additionalSuffix}"
    }
    
    scatter(vcfCompressed in allVcfCompressed) {
        File allVcfCompressedFile = vcfCompressed.vcf
    }
    Array[File]+ allVcfCompressedList = allVcfCompressedFile

    call mergeVcf.MergeCallersGetCandidate {
        input:
            pairName = pairName,
            candidateBedPrefixPath = candidateBedPrefixPath,
            allVcfCompressed = allVcfCompressed,
            allVcfCompressedList = allVcfCompressedList,
            intervalListBed = intervalListBed,
            allVcfCompressed = allVcfCompressed,
            allVcfCompressedList = allVcfCompressedList,
            minSplits = lancetSplits,
            candidateBedPaths = candidateBedPaths,
            additionalSuffix = additionalSuffix
    }
    Array[File] candidateBeds = select_all(MergeCallersGetCandidate.candidateBeds)
    scatter (candidateBed in candidateBeds) {
        String chrom = sub(basename(candidateBed), ".bed$", "")
        call calling.SliceBam as sliceBamTumor {
            input:
                chrom = chrom,
                sampleId = tumorId,
                chromBed = candidateBed,
                finalBam = tumorFinalBam,
                threads = 4,
                memoryGb = 4,
                diskSize = 10,
                serviceAccountKey = serviceAccountKey,
                gcpProject = gcpProject
        }
        
        call calling.SliceBam as sliceBamNormal {
            input:
                chrom = chrom,
                sampleId = normalId,
                chromBed = candidateBed,
                finalBam = normalFinalBam,
                threads = 4,
                memoryGb = 4,
                diskSize = 10,
                serviceAccountKey = serviceAccountKey,
                gcpProject = gcpProject
        }
        Int diskSizeLancet = (ceil( size(sliceBamNormal.lancetBamSlice.bam, "GB") + size(sliceBamTumor.lancetBamSlice.bam, "GB")) ) + 10
        call calling.LancetWGSRegional as lancetConfirm {
            input:
                chrom = chrom,
                chromBed = candidateBed,
                referenceFa = referenceFa,
                normalFinalBam = sliceBamNormal.lancetBamSlice,
                tumorFinalBam = sliceBamTumor.lancetBamSlice,
                pairName = pairName,
                lancetChromVcfPath = "~{chrom}.vcf",
                threads = 8,
                memoryGb = 40,
                diskSize = diskSizeLancet
        }

                # lancetConfirm.lancetChromVcf
        call mergeVcf.CompressIndexVcf as lancetCompressIndexVcf {
            input:
                vcf = lancetConfirm.lancetChromVcf,
                memoryGb = 4
        }

        call mergeVcf.IntersectVcfs {
            input:
                chrom = chrom,
                pairName = pairName,
                vcfCompressedLancet = lancetCompressIndexVcf.vcfCompressedIndexed ,
                vcfCompressedCandidate = MergeCallersGetCandidate.candidateVcfCompressedIndexed

        }

        #  =================================================================
        #                   Prep supporting Lancet calls
        #  =================================================================

        call mergeVcf.PrepLancetConfirm {
            input:
                pairName = pairName,
                tumorId = tumorId,
                normalId = normalId,
                tool = "lancet",
                referenceFa = referenceFa,
                callerVcf = IntersectVcfs.vcfConfirmedCandidate,
                tool = "lancet"
        }
        
        call mergeVcf.RemoveContig {
            input:
                mnvVcfPath = sub(basename(PrepLancetConfirm.mnvVcf), ".vcf$", ".removed_contig.vcf"),
                removeChromVcf = PrepLancetConfirm.mnvVcf
        }
    }
    call mergeVcf.Gatk4MergeSortVcf {
    input:
        tempVcfs = RemoveContig.removeContigVcf,
        sortedVcfPath = "~{pairName}.removed_contig.vcf.gz",
        referenceFa = referenceFa,
        gzipped = true,
        memoryGb = 16,
        diskSize = 60
    }
    call mergeVcf.MergeCallersAll as lancetMergeCallers {
        input:
            pairName = pairName,
            allVcfCompressed = [MergeCallersGetCandidate.mergedVcfCompressedIndexed, Gatk4MergeSortVcf.sortedVcf],
            allVcfCompressedList = [MergeCallersGetCandidate.mergedVcfCompressedIndexed.vcf, Gatk4MergeSortVcf.sortedVcf.vcf]
    }
    # gather split VCF filenames to avoid glob that can fail on prem
    Array[String] vcfSuffixes = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29", "30", "31", "32", "33", "34", "35", "36", "37", "38", "39", "40", "41", "42", "43", "44", "45", "46", "47", "48", "49", "50"]
    Int vcfMaxSplits = 30
    String vcfPrefix = "merged.plus.lancet"
    String vcfAdditionalSuffix = ".vcf"
    Array[String] mnvVcfPaths = ["~{vcfPrefix}.mnv~{vcfAdditionalSuffix}"]
    scatter (index in range(length(vcfSuffixes))) {
        String vcfSuffix = vcfSuffixes[index]
        String splitChromVcfPaths = "~{vcfPrefix}.~{vcfSuffix}~{vcfAdditionalSuffix}"
    }
    Array[String] splitVcfPaths = flatten([mnvVcfPaths, splitChromVcfPaths])
    call utils.SplitVcf {
        input:
            vcf=lancetMergeCallers.mergedChromVcf,
            prefix=vcfPrefix,
            diskSize = (ceil(size(lancetMergeCallers.mergedChromVcf, "GB")) * 3) + 10,
            maxRows = 1000,
            maxSplits = vcfMaxSplits,
            minSplits = 2,
            splitVcfPaths = splitVcfPaths
    }
    Array[File] splitVcfs = select_all(SplitVcf.splitVcfs)
    #  =================================================================
    #                     Merge columns
    #  =================================================================
    scatter(vcf in splitVcfs) {
        String splitId = sub(basename(vcf), ".vcf$", "")
        Int diskSize = ceil( size(tumorFinalBam.bam, "GB") + size(normalFinalBam.bam, "GB")) + 20
        call mergeVcf.PostProcessMerged {
            input:
                chrom = splitId,
                tumorId = tumorId,
                normalId = normalId,
                pairName = pairName,
                supportedChromVcf = vcf,
                normalFinalBam = normalFinalBam,
                tumorFinalBam = tumorFinalBam,
                ponFile = ponFile,
                diskSize = diskSize,
                serviceAccountKey = serviceAccountKey,
                gcpProject = gcpProject
        }
    }


    output {
            Array[File] finalChromVcf = PostProcessMerged.finalChromVcf
        }
}
