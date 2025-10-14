version 1.0

import "utils.wdl" as utils
import "../calling/calling.wdl" as calling
import "../wdl_structs.wdl"

workflow Mutect2RateTest {
    input {
        Boolean local = false
        String library

        PairInfo pairInfo
        #   mutect2
        Array[String]+ chrXCallerIntervals
        Array[File]+ chrXCallerIntervalsBedFewNodes
        File invertedIntervalListBed
        IndexedReference referenceFa
        File mutectJsonLog
        
        Boolean highMem = false
        Int diskSize = 20
        Int mutect2RateThreshold = 1000
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
        Array[String] callerIntervalsRunCloud = chrXCallerIntervals
    }
    if (local) {
        Int localCallMemoryGb = 48
        Array[File] callerIntervalsRunLocal = chrXCallerIntervalsBedFewNodes
    }
    Int callMemoryGb = select_first([localCallMemoryGb, highCallMemoryGb, lowCallMemoryGb])
    Int filterMemoryGb = select_first([highFilterMemoryGb, lowFilterMemoryGb])
    Int filterDiskSize = select_first([highFilterDiskSize, lowFilterDiskSize])
    scatter(callerInterval in select_first([callerIntervalsRunLocal, callerIntervalsRunCloud])) {
        if (library == 'WGS') {
            if (local) {
                call calling.Mutect2WgsFewNodes {
                    input:
                        callerIntervalsBed = callerInterval,
                        tumor = pairInfo.tumorId,
                        normal = pairInfo.normalId,
                        pairName = pairInfo.analysisPairId,
                        referenceFa = referenceFa,
                        normalFinalBam = pairInfo.normalFinalBam,
                        tumorFinalBam = pairInfo.tumorFinalBam,
                        memoryGb = callMemoryGb,
                        diskSize = diskSize
                }
            }
            if (!local) {
                call calling.Mutect2Wgs {
                    input:
                        chrom = callerInterval,
                        tumor = pairInfo.tumorId,
                        normal = pairInfo.normalId,
                        pairName = pairInfo.analysisPairId,
                        referenceFa = referenceFa,
                        normalFinalBam = pairInfo.normalFinalBam,
                        tumorFinalBam = pairInfo.tumorFinalBam,
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
                        tumor = pairInfo.tumorId,
                        normal = pairInfo.normalId,
                        pairName = pairInfo.analysisPairId,
                        referenceFa = referenceFa,
                        normalFinalBam = pairInfo.normalFinalBam,
                        tumorFinalBam = pairInfo.tumorFinalBam,
                        invertedIntervalListBed = invertedIntervalListBed,
                        memoryGb = callMemoryGb,
                        diskSize = diskSize
                }
            }
            if (!local) {
                call calling.Mutect2Exome {
                    input:
                        chrom = callerInterval,
                        tumor = pairInfo.tumorId,
                        normal = pairInfo.normalId,
                        pairName = pairInfo.analysisPairId,
                        referenceFa = referenceFa,
                        normalFinalBam = pairInfo.normalFinalBam,
                        tumorFinalBam = pairInfo.tumorFinalBam,
                        invertedIntervalListBed = invertedIntervalListBed,
                        memoryGb = callMemoryGb,
                        diskSize = diskSize
                }
            }
        }
        File mutect2Log = select_first([Mutect2Wgs.mutect2Log, Mutect2Exome.mutect2Log, Mutect2WgsFewNodes.mutect2Log, Mutect2ExomeFewNodes.mutect2Log])
        File mutect2ChrXRawVcfInput = select_first([Mutect2Wgs.mutect2ChromRawVcf, Mutect2Exome.mutect2ChromRawVcf, Mutect2WgsFewNodes.mutect2ChromRawVcf, Mutect2ExomeFewNodes.mutect2ChromRawVcf])
        File mutect2ChrXRawStatsInput = select_first([Mutect2Wgs.mutect2ChromRawStats, Mutect2Exome.mutect2ChromRawStats, Mutect2WgsFewNodes.mutect2ChromRawStats, Mutect2ExomeFewNodes.mutect2ChromRawStats])
    }
    call utils.GetMutect2RateReadless {
        input:
            mutect2Logs = mutect2Log,
            mutect2RateTablePath = "~{pairInfo.analysisPairId}.mutect2.rate.csv",
            mutect2RatePath = "~{pairInfo.analysisPairId}.mutect2.rate.txt"
    }
    call utils.TestMutect2Rate {
        input:
            mutect2Rate = GetMutect2RateReadless.mutect2Rate,
            mutect2RateThreshold = mutect2RateThreshold
    }
    output {
        Array[File] mutect2ChrXRawVcfs = mutect2ChrXRawVcfInput
        Array[File] mutect2ChrXRawStats = mutect2ChrXRawStatsInput
        File mutect2chrXRateTable = GetMutect2RateReadless.mutect2RateTable
        File mutect2Rate = GetMutect2RateReadless.mutect2Rate
        File? mutect2RatePassed = TestMutect2Rate.mutect2RatePassed
    }
}
