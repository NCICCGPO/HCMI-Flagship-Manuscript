version 1.0
import "../merge_vcf/merge_vcf.wdl" as mergeVcf
import "../tasks/utils.wdl" as utils
import "calling.wdl" as calling
import "../wdl_structs.wdl"


workflow Mutect2MultiGenotype {
    # command
    #   run Mutect2 caller
    input {
        Boolean local = false
        Array[String] tumors
        String normal
        File allelesVcf
        String participantId
        IndexedReference referenceFa
        Bam normalFinalBam
        Array[Bam] tumorFinalBams
        Int diskSize = 20        
        Boolean highMem = false
    }
    
    Int lowCallMemoryGb = 20
    Int lowFilterMemoryGb = 20
    Int lowFilterDiskSize = 10
    
    if (highMem) {
        Int highCallMemoryGb = 40
        Int highFilterMemoryGb = 20
        Int highFilterDiskSize = 20
    }
    
    if (local) {
        Int localCallMemoryGb = 48
    }
    Int callMemoryGb = select_first([localCallMemoryGb, highCallMemoryGb, lowCallMemoryGb])
    Int filterMemoryGb = select_first([highFilterMemoryGb, lowFilterMemoryGb])
    Int filterDiskSize = select_first([highFilterDiskSize, lowFilterDiskSize])
    
    
    # gather split VCF filenames to avoid glob that can fail on prem
    Array[String] vcfSuffixes = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29", "30", "31", "32", "33", "34", "35", "36", "37", "38", "39", "40", "41", "42", "43", "44", "45", "46", "47", "48", "49", "50"]
    Int vcfMaxSplits = 30
    String vcfPrefix = "mutect2Multisample.genotyped"
    String vcfAdditionalSuffix = ".vcf"
    Array[String] mnvVcfPaths = ["~{participantId}~{vcfPrefix}.mnv~{vcfAdditionalSuffix}"]
    scatter (index in range(length(vcfSuffixes))) {
        String vcfSuffix = vcfSuffixes[index]
        String splitChromVcfPaths = "~{participantId}~{vcfPrefix}.~{vcfSuffix}~{vcfAdditionalSuffix}"
    }
    Array[String] splitVcfPaths = flatten([mnvVcfPaths, splitChromVcfPaths])
    call utils.SplitVcf {
        input:
            vcf=allelesVcf,
            prefix="~{participantId}~{vcfPrefix}",
            diskSize = (ceil(size(allelesVcf, "GB")) * 3) + 10,
            maxRows = 1000,
            maxSplits = vcfMaxSplits,
            minSplits = 2,
            splitVcfPaths = splitVcfPaths
    }
    Array[File] splitVcfs = select_all(SplitVcf.splitVcfs)
    
    scatter (tumorFinalBam in  tumorFinalBams) {
        File tumorFinalBamList = tumorFinalBam.bam
    }

    scatter(i in range(length(splitVcfs))) {
        call mergeVcf.CompressIndexVcf as centerCallVcf {
                input:
                    vcf = splitVcfs[i],
                    memoryGb = 8
        }
        call calling.Mutect2WgsMultiGenotype {
            input:
                alleles = centerCallVcf.vcfCompressedIndexed,
                splitId = vcfSuffixes[i],
                tumors = tumors,
                normal = normal,
                participantId = participantId,
                referenceFa = referenceFa,
                normalFinalBam = normalFinalBam,
                tumorFinalBams = tumorFinalBams,
                tumorFinalBamList = tumorFinalBamList,
                memoryGb = callMemoryGb,
                diskSize = diskSize
        }
    }
    call calling.MergeMutectStats {
        input:
            pairName = participantId,
            mutect2ChromRawStats = Mutect2WgsMultiGenotype.mutect2ChromRawStats,
            memoryGb = filterMemoryGb,
            diskSize = filterDiskSize
    }
    
    call calling.LearnReadOrientationModel as learnReadOrientationModel {
        input: 
            memoryGb = 16,
            diskSize = 10,
            pairName = participantId,
            f1r2Gz =  Mutect2WgsMultiGenotype.f1r2Gz
    
    }
    
    # unfiltered
    call calling.Gatk4MergeSortVcf as unfilteredGatk4MergeSortVcf {
        input:
            sortedVcfPath = "~{participantId}.mutect2.unfiltered.sorted.vcf",
            tempChromVcfs = Mutect2WgsMultiGenotype.mutect2ChromRawVcf,
            referenceFa = referenceFa,
            memoryGb = 8,
            diskSize = 10
    }

    call calling.Mutect2FilterFfpe {
        input:
            pairName = participantId,
            referenceFa = referenceFa,
            mutect2ChromRawVcf = unfilteredGatk4MergeSortVcf.sortedVcf.vcf,
            mutect2ChromRawStats = MergeMutectStats.mergedMutect2ChromRawStats,
            outArtifacts = learnReadOrientationModel.outArtifacts,
            memoryGb = filterMemoryGb,
            diskSize = filterDiskSize
    }

    output {
        File mutect2 = Mutect2FilterFfpe.mutect2ChromVcf
        File mutect2Unfiltered = unfilteredGatk4MergeSortVcf.sortedVcf.vcf
    }
}
