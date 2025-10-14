version 1.0

import "./merge_bams.wdl" as mergeBams
import "./qc.wdl" as qc
import "../tasks/utils.wdl" as utils
import "../wdl_structs.wdl"

workflow MergeBams {
    # command
    #     merge lane level BAMs
    input {
        Boolean external = false
        #    command merge flowcell
        Array[File] laneFixmateBams
        String sampleId
        String sampleAnalysisId
        IndexedVcf MillsAnd1000G
        IndexedVcf Indels
        IndexedVcf dbsnp
        File callRegions
        #IndexedTable callRegions
        IndexedReference referenceFa
        File randomIntervals
        String qcDir = '.'       # To pass to the two QC tasks and Novosort
        # resources
        File? novosortLicense
        Int novosortMem
        Int threads
        Int additionalDiskSize = 20
        Int printReadsPreemptible = 3
    }

    Int laneFixmateBamsSize = ceil(size(laneFixmateBams, "GB"))

    if (!external) {
        if (!defined(novosortLicense)) {
            call utils.CreateBlankFile as createNovosortLicense {
                input:
                    fileId = "~{sampleId}_nonNovosort.lic_"
            }
        }

        File submittedLicense = select_first([novosortLicense, createNovosortLicense.blankFile])

        call mergeBams.NovosortMarkDup as novosort {
            input:
                laneBams = laneFixmateBams,
                sampleId = sampleId,
                logDir = qcDir,
                novosortLicense = submittedLicense,
                memoryGb = novosortMem,
                threads = threads,
                # novosort uses a lot of memory and a lot of disk.
                diskSize = (3 * laneFixmateBamsSize)
        }
    }

    if (external) {
        call mergeBams.NovosortMarkDupExternal as novosortExternal {
            input:
                laneBams = laneFixmateBams,
                sampleId = sampleId,
                memoryGb = novosortMem,
                threads = 1,
                # novosort uses a lot of memory and a lot of disk.
                diskSize = (3 * laneFixmateBamsSize)
        }
    }

    Bam mergedDedupBam = select_first([novosort.mergedDedupBam, novosortExternal.mergedDedupBam])

    # This task runs in parallel with Bqsr38. We are missing coverage check
    # tasks. The idea is that we check coverage and if it's lower than expected
    # coverage, the check task fails thereby stopping the workflow.
    String collectWgsMetricsPath = "~{qcDir}/~{sampleAnalysisId}.CollectWgsMetrics.dedup.txt"
    call qc.CollectWgsMetrics {
        input:
            inputBam = mergedDedupBam,
            sampleId = sampleId,
            outputDir = qcDir,
            collectWgsMetricsPath = collectWgsMetricsPath,
            referenceFa = referenceFa,
            randomIntervals = randomIntervals,
            diskSize = ceil(size(mergedDedupBam.bam, "GB") * 1.5) + additionalDiskSize
    }
    String  MultipleMetricsBasePreBqsrBasename = "~{qcDir}/~{sampleAnalysisId}.MultipleMetrics.dedup"
    call qc.MultipleMetricsPreBqsr {
        input:
            referenceFa = referenceFa,
            mergedDedupBam = mergedDedupBam,
            outputDir = qcDir,
             MultipleMetricsBasePreBqsrBasename =  MultipleMetricsBasePreBqsrBasename,
            sampleId = sampleId,
            diskSize = ceil(size(mergedDedupBam.bam, "GB") * 1.5) + additionalDiskSize
    }

    call mergeBams.Downsample {
        input:
            mergedDedupBam = mergedDedupBam,
            sampleId = sampleId,
            diskSize = ceil(size(mergedDedupBam.bam, "GB") * 1.5) + additionalDiskSize
    }
    call mergeBams.Bqsr38 {
        input:
            mergedDedupBam = Downsample.downsampleMergedDedupBam,
            MillsAnd1000G = MillsAnd1000G,
            referenceFa = referenceFa,
            Indels = Indels,
            dbsnp = dbsnp,
            callRegions = callRegions,
            sampleId = sampleId,
            diskSize = ceil(size(Downsample.downsampleMergedDedupBam.bam, "GB") * 1.5) + additionalDiskSize
    }

    String finalBamPath = "~{sampleAnalysisId}.final.bam"
    call mergeBams.PrintReads {
        input:
            referenceFa = referenceFa,
            mergedDedupBam = mergedDedupBam,
            finalBamPath = finalBamPath,
            recalGrp = Bqsr38.recalGrp,
            sampleId = sampleId,
            diskSize = ceil(3 * laneFixmateBamsSize)  + additionalDiskSize,
            preemptible = printReadsPreemptible
    }

    output {
        Bam finalBam = PrintReads.finalBam
        File collectWgsMetricsPreBqsr = CollectWgsMetrics.collectWgsMetrics
        File qualityDistributionPdfPreBqsr = MultipleMetricsPreBqsr.qualityDistributionPdfPreBqsr
        File qualityByCycleMetricsPreBqsr = MultipleMetricsPreBqsr.qualityByCycleMetricsPreBqsr
        File qualityByCyclePdfPreBqsr = MultipleMetricsPreBqsr.qualityByCyclePdfPreBqsr
        File qualityDistributionMetricsPreBqsr = MultipleMetricsPreBqsr.qualityDistributionMetricsPreBqsr
        File dedupLog = select_first([novosort.dedupLog, novosortExternal.dedupLog])
  }
}
