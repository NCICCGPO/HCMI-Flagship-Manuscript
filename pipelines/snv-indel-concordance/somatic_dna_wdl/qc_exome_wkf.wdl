version 1.0

import "wdl_structs.wdl"
import "pre_process/qc_exome_wkf.wdl" as qc

# breaking full pipeline to just handle the preprocessing steps to test aligners
# input is a list of samples and there fastq information, no pairing info needed

# ================== COPYRIGHT ================================================
# New York Genome Center
# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright (2021) by the New York
# Genome Center. All rights are reserved. This software is supplied without
# any warranty or guaranteed support whatsoever. The New York Genome Center
# cannot be responsible for its use, misuse, or functionality.
#
#    Jennifer M Shelton (jshelton@nygenome.org)
#    James Roche (jroche@nygenome.org)
#    Nico Robine (nrobine@nygenome.org)
#    Timothy Chu (tchu@nygenome.org)
#    Will Hooper (whooper@nygenome.org)
#    Minita Shah
#
# ================== /COPYRIGHT ===============================================


workflow QcExomeWrapper {
    input {
        IndexedReference referenceFa
        IndexedVcf gnomadBiallelic
        File chromLengths
        File intervalList
        File refSeqGeneList
        Array[SampleBamInfo] bamInfos
        File markerBedFile

        Int preprocessAdditionalDiskSize = 20
    }
    scatter (bamInfo in bamInfos) {
        call qc.QcMetricsExome {
            input:
                finalBam = bamInfo.finalBam,
                referenceFa = referenceFa,
                sampleId = bamInfo.sampleId,
                intervalList = intervalList,
                refSeqGeneList = refSeqGeneList,
                chromLengths = chromLengths,
                gnomadBiallelic = gnomadBiallelic,
                markerBedFile = markerBedFile,
                outputDir = "Sample_~{bamInfo.sampleId}/qc"
        }
    }

    output {
        # QC
        Array[File] alignmentSummaryMetrics=QcMetricsExome.alignmentSummaryMetrics
        Array[File] qualityByCyclePdf=QcMetricsExome.qualityByCyclePdf
        Array[File] baseDistributionByCycleMetrics=QcMetricsExome.baseDistributionByCycleMetrics
        Array[File] qualityByCycleMetrics=QcMetricsExome.qualityByCycleMetrics
        Array[File] baseDistributionByCyclePdf=QcMetricsExome.baseDistributionByCyclePdf
        Array[File] qualityDistributionPdf=QcMetricsExome.qualityDistributionPdf
        Array[File] qualityDistributionMetrics=QcMetricsExome.qualityDistributionMetrics
        Array[File] insertSizeHistogramPdf=QcMetricsExome.insertSizeHistogramPdf
        Array[File] insertSizeMetrics=QcMetricsExome.insertSizeMetrics
        Array[File] gcBiasMetrics=QcMetricsExome.gcBiasMetrics
        Array[File] gcBiasSummary=QcMetricsExome.gcBiasSummary
        Array[File] gcBiasPdf=QcMetricsExome.gcBiasPdf
        Array[File] flagStat=QcMetricsExome.flagStat
        Array[File] hsMetrics=QcMetricsExome.hsMetrics
        Array[File] hsMetricsPerTargetCoverage=QcMetricsExome.hsMetricsPerTargetCoverage
        Array[File] collectOxoGMetrics=QcMetricsExome.collectOxoGMetrics
        Array[File] depthOfCoverageSampleSummary=QcMetricsExome.depthOfCoverageSampleSummary
        Array[File] depthOfCoverageSampleStatistics=QcMetricsExome.depthOfCoverageSampleStatistics
        Array[File] depthOfCoverageSampleGeneSummary=QcMetricsExome.depthOfCoverageSampleGeneSummary
        Array[File] pileupsConpair=QcMetricsExome.pileupsConpair
        Array[File] contaminationTable=QcMetricsExome.contaminationTable
        Array[File] pileupsTable = QcMetricsExome.pileupsTable
    }

}
