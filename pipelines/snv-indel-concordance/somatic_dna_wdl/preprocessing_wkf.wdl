version 1.0

import "./wdl_structs.wdl"
import "tasks/utils.wdl" as utils
import "nygc_pipe_segments/nygc_germline_cancer_wkf.wdl" as germlineCancer
import "pre_process/pre_process_wkf.wdl" as preProcess
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


workflow PreprocessWrapper {
    input {
        Boolean external = false
        Boolean local = false
        Boolean sampleIsNormal = false
        Float expectedCoverage = 80
        Boolean bypassQcCheck = false
        String library = "WGS"
        Boolean production = true
        Boolean skipNormalOnlySteps = false

        BwaMem2Reference bwamem2Reference
        IndexedReference referenceFa
        File adaptersFa
        IndexedVcf MillsAnd1000G
        IndexedVcf Indels
        IndexedVcf dbsnp
        IndexedVcf gnomadBiallelic
        File bqsrCallRegions
        File chromLengths
        File hsMetricsIntervals
        File randomIntervals
        SampleInfo sampleInfo
        File markerBedFile
        File? novosortLicense
        
        File intervalListBed
        # kourami
        BwaReference kouramiReference
        File kouramiFastaGem1Index
        #fastNgsAdmix
        File fastNgsAdmixChroms
        File fastNgsAdmixContinentalSites
        File fastNgsAdmixContinentalSitesBin
        File fastNgsAdmixContinentalSitesIdx
        File fastNgsAdmixContinentalRef
        File fastNgsAdmixContinentalNind
        File fastNgsAdmixPopulationSites
        File fastNgsAdmixPopulationSitesBin
        File fastNgsAdmixPopulationSitesIdx
        File fastNgsAdmixPopulationRef
        File fastNgsAdmixPopulationNind
        # post annotation
        File cosmicCensus
        File ensemblEntrez
        # germline
        File excludeIntervalList
        Array[File] scatterIntervalsHcs
        # annotation:
        String vepGenomeBuild
        Map[String, IndexedVcf] cosmicCodingChrom
        Map[String, IndexedVcf] cosmicNoncodingChrom
        Array[String]+ listOfChroms
        Array[String]+ listOfChromsFull
        Map[String, File] vepCacheChrom
        Map[String, File] annotationsChrom
        File plugins
        File dbNSFPReplacementLogic
        File MaxEntScan
        Map[String, IndexedTable] dbNSFP4Chrom
        Map[String, IndexedTable] dbscSNVChrom
        Map[String, IndexedVcf] gnomadGenomes
        Map[String, IndexedVcf] gnomadExomes
        IndexedReference vepFastaReference
        File cancerResistanceMutations
        IndexedVcf MillsAnd1000G
        IndexedVcf omni
        IndexedVcf hapmap
        IndexedVcf onekG
        IndexedVcf dbsnp
        IndexedVcf whitelist
        IndexedVcf nygcAf
        IndexedVcf pgx
        IndexedTable rwgsPgxBed
        IndexedVcf deepIntronicsVcf
        IndexedVcf clinvarIntronicsVcf
        IndexedVcf chdWhitelistVcf

        Boolean trim = true
        Int preprocessAdditionalDiskSize = 20
        Boolean germlineHighMem = false
        # need with localization_only
        File? serviceAccountKey
        String? gcpProject
    }
    String clientSampleId = sampleInfo['listOfFastqPairs'][0].clientSampleId
    call preProcess.Preprocess {
        input:
            external = external,
            additionalDiskSize = preprocessAdditionalDiskSize,
            listOfFastqPairs = sampleInfo.listOfFastqPairs,
            trim = trim,
            adaptersFa = adaptersFa,
            sampleId = clientSampleId,
            sampleAnalysisId = sampleInfo.sampleAnalysisId,
            bwamem2Reference = bwamem2Reference,
            referenceFa = referenceFa,
            MillsAnd1000G = MillsAnd1000G,
            gnomadBiallelic = gnomadBiallelic,
            hsMetricsIntervals = hsMetricsIntervals,
            callRegions = bqsrCallRegions,
            randomIntervals = randomIntervals,
            Indels = Indels,
            dbsnp = dbsnp,
            chromLengths = chromLengths,
            markerBedFile = markerBedFile,
            novosortLicense = novosortLicense,
            local = local
    }

    # -> check normal coverage for germline pipeline
    if (!bypassQcCheck ) {
        call utils.BamQcCheck {
            input:
                wgsMetricsFile = Preprocess.collectWgsMetrics,
                expectedCoverage = expectedCoverage
        }
    }
    Boolean qcPass = select_first([BamQcCheck.coveragePass, bypassQcCheck])
    if (sampleIsNormal && qcPass && !skipNormalOnlySteps) {
        SampleBamInfo normalSampleBamInfo = object {
            #sampleId : sampleInfo.sampleAnalysisId,
            sampleId : sampleInfo.listOfFastqPairs[0].clientSampleId,
            finalBam : Preprocess.finalBam
        }
        
        call germlineCancer.NygcCancerGermline {
            input:
                production = production,
                local = local,
                normalSampleBamInfo = normalSampleBamInfo,
                analysisNormalId = sampleInfo.sampleAnalysisId,
                referenceFa = referenceFa,
                intervalListBed = intervalListBed,
                kouramiReference = kouramiReference,
                kouramiFastaGem1Index = kouramiFastaGem1Index,
                fastNgsAdmixContinentalSites = fastNgsAdmixContinentalSites,
                fastNgsAdmixContinentalSitesBin = fastNgsAdmixContinentalSitesBin,
                fastNgsAdmixContinentalSitesIdx = fastNgsAdmixContinentalSitesIdx,
                fastNgsAdmixChroms = fastNgsAdmixChroms,
                fastNgsAdmixContinentalRef = fastNgsAdmixContinentalRef,
                fastNgsAdmixContinentalNind = fastNgsAdmixContinentalNind,
                fastNgsAdmixPopulationSites = fastNgsAdmixPopulationSites,
                fastNgsAdmixPopulationSitesBin = fastNgsAdmixPopulationSitesBin,
                fastNgsAdmixPopulationSitesIdx = fastNgsAdmixPopulationSitesIdx,
                fastNgsAdmixPopulationRef = fastNgsAdmixPopulationRef,
                fastNgsAdmixPopulationNind = fastNgsAdmixPopulationNind,
                listOfChroms = listOfChroms,
                MillsAnd1000G = MillsAnd1000G,
                hapmap = hapmap,
                omni = omni,
                onekG = onekG,
                dbsnp = dbsnp,
                nygcAf = nygcAf,
                excludeIntervalList = excludeIntervalList,
                scatterIntervalsHcs = scatterIntervalsHcs,
                pgx = pgx,
                rwgsPgxBed = rwgsPgxBed,
                whitelist = whitelist,
                chdWhitelistVcf = chdWhitelistVcf,
                deepIntronicsVcf = deepIntronicsVcf,
                clinvarIntronicsVcf = clinvarIntronicsVcf,
                germlineHighMem = germlineHighMem,
                vepGenomeBuild = vepGenomeBuild,
                cosmicCodingChrom = cosmicCodingChrom,
                cosmicNoncodingChrom = cosmicNoncodingChrom,
                dbNSFPReplacementLogic = dbNSFPReplacementLogic,
                MaxEntScan = MaxEntScan,
                dbNSFP4Chrom = dbNSFP4Chrom,
                dbscSNVChrom = dbscSNVChrom,
                gnomadGenomes = gnomadGenomes,
                gnomadExomes = gnomadExomes,
                # Public
                cancerResistanceMutations = cancerResistanceMutations,
                vepCacheChrom = vepCacheChrom,
                annotationsChrom = annotationsChrom,
                plugins = plugins,
                vepFastaReference = vepFastaReference,
                # Public
                deepIntronicsVcf = deepIntronicsVcf,
                clinvarIntronicsVcf = clinvarIntronicsVcf,
                # post annotation
                cosmicCensus = cosmicCensus,
                ensemblEntrez = ensemblEntrez,
                library = library,
                serviceAccountKey = serviceAccountKey,
                gcpProject = gcpProject
        }
        PreprocessingNormalOnlyOutput normalOutput = object {
            haplotypecallerFinalFiltered: NygcCancerGermline.haplotypecallerFinalFiltered,
            filteredHaplotypecallerAnnotatedVcf: NygcCancerGermline.filteredHaplotypecallerAnnotatedVcf,
            haplotypecallerAnnotatedVcf: NygcCancerGermline.haplotypecallerAnnotatedVcf,
            kouramiResult: NygcCancerGermline.kouramiResult,
            r1HlaFastq: NygcCancerGermline.r1HlaFastq,
            r2HlaFastq: NygcCancerGermline.r2HlaFastq,
            # ancestry
            beagleFileContinental: NygcCancerGermline.beagleFileContinental,
            fastNgsAdmixQoptContinental: NygcCancerGermline.fastNgsAdmixQoptContinental,
            beagleFilePopulation: NygcCancerGermline.beagleFilePopulation,
            fastNgsAdmixQoptPopulation: NygcCancerGermline.fastNgsAdmixQoptPopulation
        }
    }

    PreprocessingOutput workflowOutput = object {
        # Bams
        finalBam: Preprocess.finalBam,
        # QC
        alignmentSummaryMetrics: Preprocess.alignmentSummaryMetrics,
        qualityByCyclePdf: Preprocess.qualityByCyclePdf,
        baseDistributionByCycleMetrics: Preprocess.baseDistributionByCycleMetrics,
        qualityByCycleMetrics: Preprocess.qualityByCycleMetrics,
        baseDistributionByCyclePdf: Preprocess.baseDistributionByCyclePdf,
        qualityDistributionPdf: Preprocess.qualityDistributionPdf,
        qualityDistributionMetrics: Preprocess.qualityDistributionMetrics,
        insertSizeHistogramPdf: Preprocess.insertSizeHistogramPdf,
        insertSizeMetrics: Preprocess.insertSizeMetrics,
        gcBiasMetrics: Preprocess.gcBiasMetrics,
        gcBiasSummary: Preprocess.gcBiasSummary,
        gcBiasPdf: Preprocess.gcBiasPdf,
        flagStat: Preprocess.flagStat,
        hsMetrics: Preprocess.hsMetrics,
        hsMetricsPerTargetCoverage: Preprocess.hsMetricsPerTargetCoverage,
        hsMetricsPerTargetCoverageAutocorr: Preprocess.hsMetricsPerTargetCoverageAutocorr,
        autocorroutput1100: Preprocess.autocorroutput1100,
        collectOxoGMetrics: Preprocess.collectOxoGMetrics,
        collectWgsMetrics: Preprocess.collectWgsMetrics,
        binestCov: Preprocess.binestCov,
        binestSex: Preprocess.binestSex,
        normCoverageByChrPng: Preprocess.normCoverageByChrPng,
        # Dedup metrics
        collectWgsMetricsPreBqsr: Preprocess.collectWgsMetricsPreBqsr,
        qualityDistributionPdfPreBqsr: Preprocess.qualityDistributionPdfPreBqsr,
        qualityByCycleMetricsPreBqsr: Preprocess.qualityByCycleMetricsPreBqsr,
        qualityByCyclePdfPreBqsr: Preprocess.qualityByCyclePdfPreBqsr,
        qualityDistributionMetricsPreBqsr: Preprocess.qualityDistributionMetricsPreBqsr,
        pileupsConpair: Preprocess.pileupsConpair,
        contaminationTable: Preprocess.contaminationTable,
        pileupsTable: Preprocess.pileupsTable
    }

    output {
        Boolean bamQCPass = qcPass
        PreprocessingOutput finalOutput = workflowOutput
        PreprocessingNormalOnlyOutput? finalNormalOutput = normalOutput
    }
}
