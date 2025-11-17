version 1.0

import "calling.wdl" as calling
import "../wdl_structs.wdl"

workflow Mutect2 {
    # command
    #   run Mutect2 caller
    input {
        Boolean testing = false
        Boolean local = false
        String library
        String tumor
        String normal
        Array[File]+ callerIntervalsBedFewNodes
        Array[String]+ callerIntervals
        Array[File] mutect2ChrXRawVcfs = []
        Array[File] mutect2ChrXRawStats = []
        String pairName
        IndexedReference referenceFa
        Bam normalFinalBam
        # Exome
        File invertedIntervalListBed
        Bam tumorFinalBam
        Int diskSize = 20
        File mutectJsonLog
        File mutectJsonLogFilter
        
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
    if (!local) {
        Array[String] callerIntervalsRunCloud = callerIntervals
    }
    if (local) {
        Int localCallMemoryGb = 48
        Array[File] callerIntervalsRunLocal = callerIntervalsBedFewNodes
    }
    Int callMemoryGb = select_first([localCallMemoryGb, highCallMemoryGb, lowCallMemoryGb])
    Int filterMemoryGb = select_first([highFilterMemoryGb, lowFilterMemoryGb])
    Int filterDiskSize = select_first([highFilterDiskSize, lowFilterDiskSize])
    
    scatter(callerInterval in select_first([callerIntervalsRunLocal, callerIntervalsRunCloud])) {
        if (library == 'WGS') {
            if (local) {
             # Mutect2WgsFewNodes
                call calling.Mutect2WgsFewNodes {
                    input:
                        callerIntervalsBed = callerInterval,
                        tumor = tumor,
                        normal = normal,
                        pairName = pairName,
                        referenceFa = referenceFa,
                        normalFinalBam = normalFinalBam,
                        tumorFinalBam = tumorFinalBam,
                        memoryGb = callMemoryGb,
                        diskSize = diskSize
                }
            } 
            if (!local) {
                call calling.Mutect2Wgs {
                    input:
                        chrom = callerInterval,
                        tumor = tumor,
                        normal = normal,
                        pairName = pairName,
                        referenceFa = referenceFa,
                        normalFinalBam = normalFinalBam,
                        tumorFinalBam = tumorFinalBam,
                        memoryGb = callMemoryGb,
                        diskSize = diskSize
                }
            }
        }
        
        if (library == 'Exome') {
            if (local) {
                call calling.Mutect2ExomeFewNodes {
                    input:
                        callerIntervalsBed = callerInterval,
                        tumor = tumor,
                        normal = normal,
                        pairName = pairName,
                        referenceFa = referenceFa,
                        normalFinalBam = normalFinalBam,
                        tumorFinalBam = tumorFinalBam,
                        invertedIntervalListBed = invertedIntervalListBed,
                        memoryGb = callMemoryGb,
                        diskSize = diskSize
                }
            }
            if (!local) {
                call calling.Mutect2Exome {
                    input:
                        chrom = callerInterval,
                        tumor = tumor,
                        normal = normal,
                        pairName = pairName,
                        referenceFa = referenceFa,
                        normalFinalBam = normalFinalBam,
                        tumorFinalBam = tumorFinalBam,
                        invertedIntervalListBed = invertedIntervalListBed,
                        memoryGb = callMemoryGb,
                        diskSize = diskSize
                }
            }
        }
        File mutect2AutosomeRawVcfInput = select_first([Mutect2Wgs.mutect2ChromRawVcf, Mutect2Exome.mutect2ChromRawVcf, Mutect2WgsFewNodes.mutect2ChromRawVcf, Mutect2ExomeFewNodes.mutect2ChromRawVcf])
        File mutect2AutosomeRawStats = select_first([Mutect2Wgs.mutect2ChromRawStats, Mutect2Exome.mutect2ChromRawStats, Mutect2WgsFewNodes.mutect2ChromRawStats, Mutect2ExomeFewNodes.mutect2ChromRawStats])
    }
    
    Array[File] mutect2AllRawVcfInput = flatten([mutect2AutosomeRawVcfInput, mutect2ChrXRawVcfs])
    Array[File] mutect2AllRawStats = flatten([mutect2AutosomeRawStats, mutect2ChrXRawStats])
    call calling.MergeMutectStats {
        input:
            pairName = pairName,
            mutect2ChromRawStats = mutect2AllRawStats,
            memoryGb = filterMemoryGb,
            diskSize = filterDiskSize
    }
    
    # unfiltered
    call calling.Gatk4MergeSortVcf as unfilteredGatk4MergeSortVcf {
        input:
            sortedVcfPath = "~{pairName}.mutect2.unfiltered.sorted.vcf",
            tempChromVcfs = mutect2AllRawVcfInput,
            referenceFa = referenceFa,
            memoryGb = 8,
            diskSize = 10
    }

    call calling.Mutect2Filter {
        input:
            pairName = pairName,
            referenceFa = referenceFa,
            mutect2ChromRawVcf = unfilteredGatk4MergeSortVcf.sortedVcf.vcf,
            mutect2ChromRawStats = MergeMutectStats.mergedMutect2ChromRawStats,
            memoryGb = filterMemoryGb,
            diskSize = filterDiskSize
    }

    call calling.AddVcfCommand as filteredAddVcfCommand {
        input:
            inVcf = Mutect2Filter.mutect2ChromVcf,
            jsonLog = mutectJsonLogFilter,
            memoryGb = 4,
            diskSize = 10
    }

    call calling.ReorderVcfColumns as filteredReorderVcfColumns {
        input:
            tumor = tumor,
            normal = normal,
            rawVcf = filteredAddVcfCommand.outVcf,
            orderedVcfPath = "~{pairName}.mutect2.vcf",
            memoryGb = 4,
            diskSize = 10
    }

    output {
        File mutect2 = filteredReorderVcfColumns.orderedVcf
    }
}
