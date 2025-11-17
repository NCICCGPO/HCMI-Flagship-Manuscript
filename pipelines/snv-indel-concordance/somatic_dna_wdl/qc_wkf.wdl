version 1.0

import "wdl_structs.wdl"
import "pre_process/qc_wkf.wdl" as qc

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


workflow QcWrapper {
    input {
        IndexedReference referenceFa
        IndexedVcf gnomadBiallelic
        File chromLengths
        File hsMetricsIntervals
        File randomIntervals
        Array[SampleBamInfo] bamInfos
        File markerBedFile

        Boolean local = true
        Int preprocessAdditionalDiskSize = 20
    }
    
    scatter (bamInfo in bamInfos) {
        call qc.QcMetrics {
            input:
                finalBam = bamInfo.finalBam,
                referenceFa = referenceFa,
                sampleId = bamInfo.sampleId,
                hsMetricsIntervals = hsMetricsIntervals,
                randomIntervals = randomIntervals,
                chromLengths = chromLengths,
                gnomadBiallelic = gnomadBiallelic,
                markerBedFile = markerBedFile,
                outputDir = "Sample_~{bamInfo.sampleId}/qc",
                local = local
        }
    }

    output {
        # QC
        Array[File] alignmentSummaryMetrics=QcMetrics.alignmentSummaryMetrics
        Array[File] qualityByCyclePdf=QcMetrics.qualityByCyclePdf
        Array[File] baseDistributionByCycleMetrics=QcMetrics.baseDistributionByCycleMetrics
        Array[File] qualityByCycleMetrics=QcMetrics.qualityByCycleMetrics
        Array[File] baseDistributionByCyclePdf=QcMetrics.baseDistributionByCyclePdf
        Array[File] qualityDistributionPdf=QcMetrics.qualityDistributionPdf
        Array[File] qualityDistributionMetrics=QcMetrics.qualityDistributionMetrics
        Array[File] insertSizeHistogramPdf=QcMetrics.insertSizeHistogramPdf
        Array[File] insertSizeMetrics=QcMetrics.insertSizeMetrics
        Array[File] gcBiasMetrics=QcMetrics.gcBiasMetrics
        Array[File] gcBiasSummary=QcMetrics.gcBiasSummary
        Array[File] gcBiasPdf=QcMetrics.gcBiasPdf
        Array[File] flagStat=QcMetrics.flagStat
        Array[File] hsMetrics=QcMetrics.hsMetrics
        Array[File] hsMetricsPerTargetCoverage=QcMetrics.hsMetricsPerTargetCoverage
        Array[File] hsMetricsPerTargetCoverageAutocorr=QcMetrics.hsMetricsPerTargetCoverageAutocorr
        Array[File] autocorroutput1100=QcMetrics.autocorroutput1100
        Array[File] collectOxoGMetrics=QcMetrics.collectOxoGMetrics
        Array[File] collectWgsMetrics=QcMetrics.collectWgsMetrics
        Array[File] binestCov=QcMetrics.binestCov
        Array[File] binestSex=QcMetrics.binestSex
        Array[File] normCoverageByChrPng=QcMetrics.normCoverageByChrPng
        Array[File] pileupsConpair=QcMetrics.pileupsConpair
        Array[File] contaminationTable=QcMetrics.contaminationTable
        Array[File] pileupsTable=QcMetrics.pileupsTable
    }

}
