version 1.0


import "qc.wdl" as qc
import "../wdl_structs.wdl"

workflow QcMetricsExome {
    # command
    input {
        Bam finalBam
        IndexedReference referenceFa
        String sampleId
        File intervalList
        File refSeqGeneList
        File chromLengths
        IndexedVcf gnomadBiallelic
        File markerBedFile
        String outputDir = "."
    }

    Int additionalDiskSize = 50
    Int diskSize = ceil((size(finalBam.bam, "GB") + size(finalBam.bamIndex, "GB"))) +
                      additionalDiskSize

    call qc.MultipleMetrics {
        input:
            referenceFa = referenceFa,
            finalBam = finalBam,
            sampleId = sampleId,
            outputDir = outputDir,
            diskSize = diskSize
    }

    call qc.CollectGcBiasMetrics {
        input:
            referenceFa = referenceFa,
            finalBam = finalBam,
            sampleId = sampleId,
            outputDir = outputDir,
            diskSize = diskSize,
    }

    call qc.Flagstat {
        input:
            referenceFa = referenceFa,
            finalBam = finalBam,
            sampleId = sampleId,
            outputDir = outputDir,
            diskSize = diskSize

    }

    call qc.HsMetricsExome {
        input:
            referenceFa = referenceFa,
            intervalList = intervalList,
            finalBam = finalBam,
            sampleId = sampleId,
            outputDir = outputDir,
            diskSize = diskSize
    }

    call qc.DepthOfCoverageExome {
        input:
            referenceFa = referenceFa,
            finalBam = finalBam,
            sampleId = sampleId,
            outputDir = outputDir,
            intervalList = intervalList,
            refSeqGeneList = refSeqGeneList,
            diskSize = diskSize
    }
    
    call qc.CollectOxoGMetricsExome {
        input:
            referenceFa = referenceFa,
            intervalList = intervalList,
            finalBam = finalBam,
            sampleId = sampleId,
            outputDir = outputDir,
            diskSize = diskSize
    }

    call qc.Pileup {
        input:
            finalBam = finalBam,
            sampleId = sampleId,
            outputDir = outputDir,
            gnomadBiallelic = gnomadBiallelic,
            diskSize = diskSize
    }

    call qc.CalculateContamination {
        input:
            sampleId = sampleId,
            pileupsTable = Pileup.pileupsTable,
            outputDir = outputDir,
            diskSize = diskSize
    }

    call qc.ConpairPileup {
        input:
            markerBedFile = markerBedFile,
            referenceFa = referenceFa,
            finalBam = finalBam,
            sampleId = sampleId,
            memoryGb = 4,
            threads = 1,
            diskSize = diskSize
    }

    output {
        File alignmentSummaryMetrics = MultipleMetrics.alignmentSummaryMetrics
        File qualityByCyclePdf = MultipleMetrics.qualityByCyclePdf
        File baseDistributionByCycleMetrics = MultipleMetrics.baseDistributionByCycleMetrics
        File qualityByCycleMetrics = MultipleMetrics.qualityByCycleMetrics
        File baseDistributionByCyclePdf = MultipleMetrics.baseDistributionByCyclePdf
        File qualityDistributionPdf = MultipleMetrics.qualityDistributionPdf
        File qualityDistributionMetrics = MultipleMetrics.qualityDistributionMetrics
        File insertSizeHistogramPdf = MultipleMetrics.insertSizeHistogramPdf
        File insertSizeMetrics = MultipleMetrics.insertSizeMetrics
        File gcBiasMetrics = CollectGcBiasMetrics.gcBiasMetrics
        File gcBiasSummary = CollectGcBiasMetrics.gcBiasSummary
        File gcBiasPdf = CollectGcBiasMetrics.gcBiasPdf
        File flagStat = Flagstat.flagStat
        File hsMetrics = HsMetricsExome.hsMetrics
        File hsMetricsPerTargetCoverage = HsMetricsExome.hsMetricsPerTargetCoverage
        File collectOxoGMetrics = CollectOxoGMetricsExome.collectOxoGMetrics
        File pileupsConpair = ConpairPileup.pileupsConpair
        File depthOfCoverageSampleSummary = DepthOfCoverageExome.depthOfCoverageSampleSummary
        File depthOfCoverageSampleStatistics = DepthOfCoverageExome.depthOfCoverageSampleStatistics
        File depthOfCoverageSampleGeneSummary = DepthOfCoverageExome.depthOfCoverageSampleGeneSummary
        File contaminationTable = CalculateContamination.contaminationTable
        File pileupsTable = Pileup.pileupsTable
    }

}
