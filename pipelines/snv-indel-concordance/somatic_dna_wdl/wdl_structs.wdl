version 1.0

# ================== COPYRIGHT ================================================
# New York Genome Center
# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright (2021) by the New York
# Genome Center. All rights are reserved. This software is supplied without
# any warranty or guaranteed support whatsoever. The New York Genome Center
# cannot be responsible for its use, misuse, or functionality.
#
#    Jennifer M Shelton (jshelton@nygenome.org)
#    Nico Robine (nrobine@nygenome.org)
#    Minita Shah (mshah@nygenome.org)
#    Timothy Chu (tchu@nygenome.org)
#    Will Hooper (whooper@nygenome.org)
#
# ================== /COPYRIGHT ===============================================

struct IndexedVcf {
    File vcf
    File index
}

struct IndexedTable {
    File table
    File index
}

struct Bam {
    File bam
    File bamIndex
    String? md5sum
}

struct Cram {
    File cram
    File cramIndex
    String? md5sum
}

struct BwaReference {
    File fasta
    File sa
    File pac
    File bwt
    File ann
    File? alt
    File amb
    File? dict
    File? index
}

struct BwaMem2Reference {
    File fasta
    File bwt2bit
    File pac
    File ann
    File amb
    File num
    File? alt
    File? dict
    File? index
}

struct Minimap2Reference {
    File fasta
    File? alt
    File? dict
    File? index
}

struct IndexedReference {
    File fasta
    File? dict
    File index
    String? httpFasta
    String? httpIndex
}

struct Fastqs {
    File fastqR1
    String? md5sumR1
    File fastqR2
    String? md5sumR2
    String readGroupId              # file name prefix
    String clientSampleId           # @RG.SM tag
    String? limsLibraryName         # @RG.LB tag
    String readGroupPlatformUnit    # @RG.PU tag
}

struct SampleInfo {
    String sampleAnalysisId
    Array[Fastqs] listOfFastqPairs
    Float expectedCoverage
    Boolean? skipCoverageCheck
}

struct SampleBamInfo {
    String sampleId
    Bam finalBam
}

struct SampleCramInfo {
    String sampleId
    Cram finalCram
}

struct PreMergedPairVcfInfo {
    String pairId
    File filteredMantaSV
    File strelka2Snv
    File strelka2Indel
    File mutect2
    File lancet
    File? svabaSv
    File? svabaIndel
    String tumor
    String normal
    Bam tumorFinalBam
    Bam normalFinalBam
}

struct PairRawVcfInfo {
    String pairId
    File? mergedVcf
    File? mainVcf
    File? supplementalVcf
    File filteredMantaSV
    File strelka2Snv
    File strelka2Indel
    File mutect2
    File lancet
    File? svabaSv
    File? svabaIndel
    IndexedVcf gridssVcf
    File bicseq2Png
    File bicseq2
    String tumor
    String normal
    Bam tumorFinalBam
    Bam normalFinalBam
}

struct MergedPairVcfInfo {
    String pairId
    String tumor
    String normal
    File unannotatedVcf
}

struct PairVcfInfo {
    String pairId
    String tumor
    String normal
    File mainVcf
    File supplementalVcf
    File vcfAnnotatedTxt
    File maf
}

struct FinalVcfPairInfo {
    String pairId
    String tumor
    String normal
    File mainVcf
    File supplementalVcf
    File filteredMantaSV
    File strelka2Snv
    File strelka2Indel
    File mutect2
    File lancet
    File? svabaSv
    File? svabaIndel
    IndexedVcf gridssVcf
    File bicseq2Png
    File bicseq2
    File cnvAnnotatedFinalBed
    File cnvAnnotatedSupplementalBed
    File svFinalBedPe
    File svHighConfidenceFinalBedPe
    File svSupplementalBedPe
    File svHighConfidenceSupplementalBedPe
}

struct FinalPairInfo {
    String pairId
    String tumor
    String normal
    File mainVcf
    File supplementalVcf
    File filteredMantaSV
    File strelka2Snv
    File strelka2Indel
    File mutect2
    File lancet
    File? svabaSv
    File? svabaIndel
    IndexedVcf gridssVcf
    File bicseq2Png
    File bicseq2
    File cnvAnnotatedFinalBed
    File cnvAnnotatedSupplementalBed
    File svFinalBedPe
    File svHighConfidenceFinalBedPe
    File svSupplementalBedPe
    File svHighConfidenceSupplementalBedPe
    Bam tumorFinalBam
    Bam normalFinalBam
}

struct PairRelationship {
    String pairId
    String tumorPrefix
    String normalPrefix
    String tumorId    #ID in SM tag
    String normalId
}

struct PairInfo {
    # internal ID is analysis ID (may be different from client ID and is used for filenaming only)
    String pairId
    String analysisPairId
    Bam tumorFinalBam
    Bam normalFinalBam
    String tumorId           # ID in SM tag
    String normalId
    String analysisTumorId           # ID in SM tag
    String analysisNormalId
}

struct PairCramInfo {
    # internal ID is analysis ID (may be different from client ID and is used for filenaming only)
    String pairId
    String analysisPairId
    Cram tumorFinalCram
    Cram normalFinalCram
    String tumor
    String normal    String analysisTumorId           # ID in SM tag
    String analysisTumorId           # ID in SM tag
    String analysisNormalId
}

struct SomaticWorkflowOutput {
    # alignment and calling results (calling results may not exist if qc failed)
    # SNV INDELs CNV SV and BAM output
    Array[FinalVcfPairInfo?] finalPairInfo

    # MSI
    Array[File?] mantisWxsKmerCountsFinal
    Array[File?] mantisWxsKmerCountsFiltered
    Array[File?] mantisExomeTxt
    Array[File?] mantisStatusFinal
    # SIGs
    Array[File?] deconstructSigSigs
    Array[File?] deconstructSigCounts
    Array[File?] deconstructSigSigInput
    Array[File?] deconstructSigReconstructed
    Array[File?] deconstructSigDiff

    # ancestry
    Array[File?] beagleFileContinental
    Array[File?] fastNgsAdmixQoptContinental
    Array[File?] beagleFilePopulation
    Array[File?] fastNgsAdmixQoptPopulation

    # Preprocessing output
    Array[Cram] finalCrams

    # QC
    Array[File] alignmentSummaryMetrics
    Array[File] qualityByCyclePdf
    Array[File] baseDistributionByCycleMetrics
    Array[File] qualityByCycleMetrics
    Array[File] baseDistributionByCyclePdf
    Array[File] qualityDistributionPdf
    Array[File] qualityDistributionMetrics
    Array[File] insertSizeHistogramPdf
    Array[File] insertSizeMetrics
    Array[File] gcBiasMetrics
    Array[File] gcBiasSummary
    Array[File] gcBiasPdf
    Array[File] flagStat
    Array[File] hsMetrics
    Array[File] hsMetricsPerTargetCoverage
    Array[File] hsMetricsPerTargetCoverageAutocorr
    Array[File] autocorroutput1100
    Array[File] collectOxoGMetrics
    Array[File] collectWgsMetrics
    Array[File] binestCov
    Array[File] binestSex
    Array[File] normCoverageByChrPng
    # Dedup metrics
    Array[File] collectWgsMetricsPreBqsr
    Array[File] qualityDistributionPdfPreBqsr
    Array[File] qualityByCycleMetricsPreBqsr
    Array[File] qualityByCyclePdfPreBqsr
    Array[File] qualityDistributionMetricsPreBqsr
    Array[File] dedupLog
    Array[File] contaminationTable
    Array[File] pileupsTable

    # Conpair
    Array[File] pileupsConpair   # per sample, not per pair
    Array[File] concordanceAll
    Array[File] concordanceHomoz
    Array[File] contamination


    # Germline
    Array[File?] kouramiResult
    Array[IndexedVcf?] haplotypecallerVcf
    Array[IndexedVcf?] haplotypecallerFinalFiltered
    Array[File?] filteredHaplotypecallerAnnotatedVcf
    Array[File?] haplotypecallerAnnotatedVcf
    Array[File?] alleleCountsTxt
}

struct PreprocessingOutput {
    # Preprocessing output.
    Bam finalBam

    # QC
    File alignmentSummaryMetrics
    File qualityByCyclePdf
    File baseDistributionByCycleMetrics
    File qualityByCycleMetrics
    File baseDistributionByCyclePdf
    File qualityDistributionPdf
    File qualityDistributionMetrics
    File insertSizeHistogramPdf
    File insertSizeMetrics
    File gcBiasMetrics
    File gcBiasSummary
    File gcBiasPdf
    File flagStat
    File hsMetrics
    File hsMetricsPerTargetCoverage
    File hsMetricsPerTargetCoverageAutocorr
    File autocorroutput1100
    File collectOxoGMetrics
    File collectWgsMetrics
    File binestCov
    File binestSex
    File normCoverageByChrPng
    File contaminationTable
    File pileupsTable

    # Dedup metrics
    File collectWgsMetricsPreBqsr
    File qualityDistributionPdfPreBqsr
    File qualityByCycleMetricsPreBqsr
    File qualityByCyclePdfPreBqsr
    File qualityDistributionMetricsPreBqsr

    # Conpair
    File pileupsConpair
}

struct PreprocessingNormalOnlyOutput {
    IndexedVcf? haplotypecallerFinalFiltered
    File? filteredHaplotypecallerAnnotatedVcf
    File? haplotypecallerAnnotatedVcf
    File? kouramiResult
    File? r1HlaFastq
    File? r2HlaFastq
    # ancestry
    File? beagleFileContinental
    File? fastNgsAdmixQoptContinental
    File? beagleFilePopulation
    File? fastNgsAdmixQoptPopulation
}
