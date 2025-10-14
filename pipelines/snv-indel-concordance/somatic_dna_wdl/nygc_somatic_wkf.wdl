version 1.0
import "wdl_structs.wdl"
import "tasks/utils.wdl" as utils
import "tasks/mutect2_test_wkf.wdl" as mutect2Test
import "nygc_pipe_segments/nygc_paired_wkf.wdl" as miscPairedTasks
import "nygc_pipe_segments/nygc_germline_cancer_wkf.wdl" as germlineCancer
import "nygc_pipe_segments/nygc_post_processing_wkf.wdl" as postProcessing
import "nygc_pipe_segments/nygc_calling_wkf.wdl" as calling
import "tasks/bam_cram_conversion.wdl" as cramConversion
import "nygc_pipe_segments/nygc_baf_wkf.wdl" as baf
import "pre_process/qc.wdl" as qc

# ================== COPYRIGHT ================================================
# New York Genome Center
# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright (2024) by the New York
# Genome Center. All rights are reserved. This software is supplied without
# any warranty or guaranteed support whatsoever. The New York Genome Center
# cannot be responsible for its use, misuse, or functionality.
#
#    Jennifer M Shelton (jshelton@nygenome.org)
#    James Roche (jroche@nygenome.org)
#    Nico Robine (nrobine@nygenome.org)
#    Timothy Chu (tchu@nygenome.org)
#    Will Hooper (whooper@nygenome.org)
#
# ================== /COPYRIGHT ===============================================
        
workflow NygcSomaticPipeline {
    input {
        Boolean tumorSkipCoverageCheck = false
        Boolean normalSkipCoverageCheck = false
        Boolean skipNormalOnlySteps = false
        Boolean trim = true
        Boolean production = true
        Boolean external = false
        Boolean local = false
        String fileType = "bam"
        Boolean bypassQcCheck = false
        Boolean assumeCallingSpeedy = false
        Int tumorExpectedCoverage = 80
        Int normalExpectedCoverage = 40
        String library = "WGS"
        # kourami
        BwaReference kouramiReference
        File kouramiFastaGem1Index
        File teloTargetIndexGem
        # mantis
        File mantisBed
        File intervalListBed
        IndexedReference referenceFa
        File markerTxtFile
        File markerBedFile
        Float maxMismatches = 0.2
        File adaptersFa
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
        Array[File] allCallerIntervalsBedFewNodes
        # annotation:
        String vepGenomeBuild
        Map[String, IndexedVcf] cosmicCodingChrom
        Map[String, IndexedVcf] cosmicNoncodingChrom
        Array[String]+ listOfChroms
        Array[String]+ listOfChromsFull
        Array[String]+ callerIntervals
        Array[File]+ callerIntervalsBedFewNodes
        Array[String]+ exomeCallerIntervals = callerIntervals
        Array[String]+ chrXCallerIntervals
        Array[File]+ chrXCallerIntervalsBedFewNodes
        Array[String]+ exomeChrXCallerIntervals = chrXCallerIntervals
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
        File invertedIntervalListBed
        #   Manta
        IndexedTable callRegions
        #   Lancet
        Array[File] splitBedsWgs
        Array[File] splitBedsWgsFewNodes
        Map[String, File] chromBeds
        #   BicSeq2
        Int readLength
        Int coordReadLength
        Map[Int, Map[String, File]] uniqCoords
        File bicseq2ConfigFile
        File bicseq2SegConfigFile
        Map[String, File] chromFastas
        Int tumorMedianInsertSize = 0
        Int normalMedianInsertSize = 0
        # Gridss
        BwaReference bwaReference
        String bsGenome
        File ponTarGz
        Array[File] gridssAdditionalReference
        # Strelka2
        File configureStrelkaSomaticWorkflow
        File lancetJsonLog
        File mantaJsonLog
        File strelkaJsonLog
        File mutectJsonLog
        File mutectJsonLogFilter
        File ponWGSFile
        File ponExomeFile
        IndexedVcf germFile
        # annotate cnv
        File cytoBand
        File dgv
        File thousandG
        File cosmicUniqueBed
        File cancerCensusBed
        File ensemblUniqueBed
        # gap,DGV,1000G,PON,COSMIC
        File gap
        File dgvBedpe
        File thousandGVcf
        File svPon
        File cosmicBedPe
        # signatures
        File cosmicSigs
        File cosmicSbsReference
        File cosmicDbsReference
        File cosmicIdReference

        # -> Patient-specific inputs
        PairInfo initialPairInfo
        Int purity = 1
        File? tumorCollectWgsMetrics
        File? normalCollectWgsMetrics
        
        # -> Preexisting files
        Array[File?] mutect2ChrXRawVcfsPreexist = []
        Array[File?] mutect2ChrXRawStatsPreexist = []
        File? normalAlignmentSummaryMetricsPreexist
        File? normalAlleleCountsTxtPreexist
        File? normalBaseDistributionByCycleMetricsPreexist
        File? normalBaseDistributionByCyclePdfPreexist
        File? normalQualityByCycleMetricsPreexist
        File? normalQualityByCyclePdfPreexist
        File? normalQualityDistributionMetricsPreexist
        File? normalQualityDistributionPdfPreexist
        File? normalInsertSizeHistogramPdfPreexist
        File? normalInsertSizeMetricsPreexist

        File? tumorAlignmentSummaryMetricsPreexist
        File? tumorAlleleCountsTxtPreexist
        File? tumorBaseDistributionByCycleMetricsPreexist
        File? tumorBaseDistributionByCyclePdfPreexist
        File? tumorQualityByCycleMetricsPreexist
        File? tumorQualityByCyclePdfPreexist
        File? tumorQualityDistributionMetricsPreexist
        File? tumorQualityDistributionPdfPreexist
        File? tumorInsertSizeHistogramPdfPreexist
        File? tumorInsertSizeMetricsPreexist

        File? alleleCountsTxtPreexist
        File? conpairPileupNormalPreexist
        File? conpairPileupTumorPreexist
        IndexedVcf? haplotypecallerFinalFilteredPreexist
        IndexedVcf? haplotypecallerVcfPreexist
        File? filteredHaplotypecallerAnnotatedVcfPreexist
        File? haplotypecallerAnnotatedVcfPreexist
        File? kouramiResultPreexist
        File? r1HlaFastqPreexist
        File? r2HlaFastqPreexist
        File? beagleFileContinentalPreexist
        File? fastNgsAdmixQoptContinentalPreexist
        File? beagleFilePopulationPreexist
        File? fastNgsAdmixQoptPopulationPreexist
        File? mutect2Preexist
        # Manta
        IndexedVcf? candidateSmallIndelsPreexist
        IndexedVcf? diploidSVPreexist
        IndexedVcf?  somaticSVPreexist
        IndexedVcf?  candidateSVPreexist
        File? unfilteredMantaSVPreexist
        File? filteredMantaSVPreexist
        # Strelka2
        IndexedVcf? strelka2SnvsPreexist
        IndexedVcf? strelka2IndelsPreexist
        File? strelka2SnvPreexist
        File? strelka2IndelPreexist
        # Lancet
        File? lancetPreexist
        # Gridss
        IndexedVcf? gridssVcfPreexist
        File? gridssUnfilteredVcfChromsPreexist
        # Bicseq2
        File? bicseq2PngPreexist
        File? bicseq2Preexist
        File? mergedVcfPreexist
        File? mainVcfPreexist
        File? supplementalVcfPreexist
        File? vcfAnnotatedTxtPreexist
        File? mafPreexist
        File? cnvAnnotatedFinalBedPreexist
        File? cnvAnnotatedSupplementalBedPreexist
        File? svFinalBedPePreexist
        File? svHighConfidenceFinalBedPePreexist
        File? svSupplementalBedPePreexist
        File? svHighConfidenceSupplementalBedPePreexist
        # Deconstructsigs
        File? sigsPreexist
        File? countsPreexist
        File? sig_inputPreexist
        File? reconstructedPreexist
        File? diffPreexist
        # Musical
        File? filteredVcfPreexist
        File? outputSbsMatrixPreexist
        File? outputDbsMatrixPreexist
        File? outputIdMatrixPreexist
        File? sigsSbsPreexist
        File? sigsDbsPreexist
        File? sigsIdPreexist
        # Gridss resources need a lot of fine grained control
        Int gridssPreMemoryGb = 60
        Int gridssFilterMemoryGb = 32
        Int annotateBicSeq2CnvMem = 36
        Boolean gridssHighMem = false
        # bicseq2 option
        Int lambda = 4
        Boolean mantaHighMem = false
        Boolean mutect2HighMem = false
        Boolean germlineHighMem = false
        
        # need for docker with gcloud
        File? serviceAccountKey
        String? gcpProject
    }
    if (local) {
        # allCallerIntervalsBedFewNodes
        Boolean LocalAssumeCallingSpeedy = true
        Array[File] scatterIntervalsHcsFewNodes = allCallerIntervalsBedFewNodes
    }
    if (!local) {
        Array[File] scatterIntervalsHcsCloud = scatterIntervalsHcs
    }
    Array[File] scatterIntervalsHcsRun = select_first([scatterIntervalsHcsFewNodes, scatterIntervalsHcsCloud])
    Boolean runAssumeCallingSpeedy = select_first([LocalAssumeCallingSpeedy, assumeCallingSpeedy])
    if (library == "Exome") {
        Array[String]+ exomeRunCallerIntervals = exomeCallerIntervals
        Array[String]+ exomeChrXRunCallerIntervals = exomeChrXCallerIntervals
    }
    if (library == "WGS") {
        Array[String]+ wgsRunCallerIntervals = callerIntervals
        Array[String]+ wgsChrXRunCallerIntervals = chrXCallerIntervals
    }
    Array[String]+ runCallerIntervals = select_first([exomeRunCallerIntervals, wgsRunCallerIntervals])
    Array[String]+ chrXRunCallerIntervals = select_first([exomeChrXRunCallerIntervals, wgsChrXRunCallerIntervals])
    # -> convert Cram to bam if needed
    if (fileType == "cram") {
        Cram tumorCram = object {
                cram : initialPairInfo.tumorFinalBam.bam,
                cramIndex : initialPairInfo.tumorFinalBam.bai
            }
        call cramConversion.SamtoolsCramToBam as tumorCramToBam {
            input:
                inputCram = tumorCram,
                referenceFa = referenceFa,
                sampleId = initialPairInfo.tumorId,
                diskSize = (ceil(size(initialPairInfo.tumorFinalBam.bam, "GB") * 3)) + 20
        }
        Cram normalCram = object {
                cram : initialPairInfo.normalFinalBam.bam,
                cramIndex : initialPairInfo.normalFinalBam.bai
            }
        call cramConversion.SamtoolsCramToBam as normalCramToBam {
            input:
                inputCram = normalCram,
                referenceFa = referenceFa,
                sampleId = initialPairInfo.normalId,
                diskSize = (ceil(size(initialPairInfo.normalFinalBam.bam, "GB") * 3)) + 20
        }
        Bam tumorFinalBamConverted = tumorCramToBam.bamInfo.finalBam
        Bam normalFinalBamConverted = normalCramToBam.bamInfo.finalBam
    }
    Bam tumorFinalBam = select_first([tumorFinalBamConverted, initialPairInfo.tumorFinalBam])
    Bam normalFinalBam = select_first([normalFinalBamConverted, initialPairInfo.normalFinalBam])
    PairInfo pairInfo = object {
            pairId : initialPairInfo.pairId,
            analysisPairId : initialPairInfo.analysisPairId,
            tumorFinalBam : tumorFinalBam, 
            normalFinalBam : normalFinalBam,
            tumorId : initialPairInfo.tumorId,
            normalId : initialPairInfo.normalId,
            analysisTumorId : initialPairInfo.analysisTumorId,
            analysisNormalId : initialPairInfo.analysisNormalId
    }
    SampleBamInfo normalSampleBamInfo = object {
        sampleId : initialPairInfo.normalId,
        finalBam : normalFinalBam
    }
    # -> check normal coverage for germline pipeline
    if (!bypassQcCheck && !normalSkipCoverageCheck && defined(normalCollectWgsMetrics)) {
        call utils.BamQcCheck as normalBamQcCheck {
            input:
                wgsMetricsFile = select_first([normalCollectWgsMetrics]),
                expectedCoverage = normalExpectedCoverage
        }
    }
    # confirm BAM QC
    Boolean normalCoveragePass = select_first([normalBamQcCheck.coveragePass, false])
    if (normalCoveragePass && costPass) {
        Boolean successBamQCGood = true
    }
    Boolean successBamQC = select_first([successBamQCGood, false])
    if (bypassQcCheck || normalSkipCoverageCheck || successBamQC ) {
        Boolean bamQCPass = true
    }
    Boolean bamQCRun = select_first([bamQCPass, false])
    # will run short tasks in NygcPairedTasks and then test Qc
    if (!bypassQcCheck && defined(tumorCollectWgsMetrics) && defined(normalCollectWgsMetrics)) {
        call utils.SomaticQcCheck {
            input:
                pairId = pairInfo.pairId,
                tumorWgsMetricsFile = select_first([tumorCollectWgsMetrics]),
                tumorExpectedCoverage = tumorExpectedCoverage,
                tumorSkipCoverageCheck = tumorSkipCoverageCheck,
                normalWgsMetricsFile = select_first([normalCollectWgsMetrics]),
                normalExpectedCoverage = normalExpectedCoverage,
                normalSkipCoverageCheck = normalSkipCoverageCheck,
                concordanceFile = NygcPairedTasks.concordanceHomoz,
                contaminationFile = NygcPairedTasks.contamination
        }
    }
    Boolean somaticQcPass = select_first([SomaticQcCheck.qcPass, bypassQcCheck])
    if (!defined(mutect2Preexist)) {
        if (length(mutect2ChrXRawVcfsPreexist) == 0) {
            call mutect2Test.Mutect2RateTest {
                input:
                    local = local,
                    library = library,
                    pairInfo = pairInfo,
                    chrXCallerIntervalsBedFewNodes = chrXCallerIntervalsBedFewNodes,
                    chrXCallerIntervals = chrXRunCallerIntervals,
                    invertedIntervalListBed = invertedIntervalListBed,
                    referenceFa = referenceFa,
                    mutectJsonLog = mutectJsonLog,
                    highMem = mutect2HighMem
            }
            if (!runAssumeCallingSpeedy) {
                call utils.TestMutect2Rate {
                    input:
                        mutect2Rate = Mutect2RateTest.mutect2Rate,
                        mutect2RateThreshold = 1000
                }
                if (defined(TestMutect2Rate.mutect2RatePassed)) {
                    Boolean fastCostPass = true
                }
            }
        }
    }
    Boolean costPassResult = select_first([fastCostPass, false])
    if (runAssumeCallingSpeedy || costPassResult ) {
        Boolean successCostPass = true
    }
    Boolean costPass = select_first([successCostPass, false])
    if (somaticQcPass && costPass) {
        Boolean costAndQcCheck = true
    }
    Boolean costAndQcPass = select_first([costAndQcCheck, bypassQcCheck])
    Array[File] mutect2ChrXRawVcfs = select_first([Mutect2RateTest.mutect2ChrXRawVcfs, select_all(mutect2ChrXRawVcfsPreexist)])
    Array[File] mutect2ChrXRawStats = select_first([Mutect2RateTest.mutect2ChrXRawStats, select_all(mutect2ChrXRawStatsPreexist)])
    # -> get insert size estimate if not provided
    Int tumorDiskSize = ceil(size(pairInfo.tumorFinalBam.bam, "GB")) + 30
    Int normalDiskSize = ceil(size(pairInfo.normalFinalBam.bam, "GB")) + 30
    if (normalMedianInsertSize == 0) {
        if (!defined(tumorInsertSizeMetricsPreexist)) {
            # tumor insert size
            String tumorMultipleMetricsBase = "./~{pairInfo.analysisTumorId}.MultipleMetrics"
            call qc.MultipleMetrics as tumorMultipleMetrics {
                input:
                    referenceFa = referenceFa,
                    finalBam = pairInfo.tumorFinalBam,
                    sampleId = pairInfo.tumorId,
                    MultipleMetricsBase = tumorMultipleMetricsBase,
                    diskSize = tumorDiskSize
            }
        }
        if (!defined(normalInsertSizeMetricsPreexist)) {
            # normal insert size
            String normalMultipleMetricsBase = "./~{pairInfo.analysisNormalId}.MultipleMetrics"
            call qc.MultipleMetrics as normalMultipleMetrics {
                input:
                    referenceFa = referenceFa,
                    finalBam = pairInfo.normalFinalBam,
                    sampleId = pairInfo.normalId,
                    MultipleMetricsBase = normalMultipleMetricsBase,
                    diskSize = normalDiskSize
            }
        }
    }
    if (bamQCRun) {
        if (!skipNormalOnlySteps) {
            call germlineCancer.NygcCancerGermline {
                input:
                    production = production,
                    local = local,
                    analysisNormalId = pairInfo.analysisNormalId,
                    normalSampleBamInfo = normalSampleBamInfo,
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
                    scatterIntervalsHcs = scatterIntervalsHcsRun,
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
        }
    }
    call miscPairedTasks.NygcPairedTasks {
        input:
            referenceFa = referenceFa,
            mantisBed = mantisBed,
            intervalListBed = intervalListBed,
            markerBedFile = markerBedFile,
            markerTxtFile = markerTxtFile,
            pairInfo = pairInfo,
            trim = trim,
            teloTargetIndexGem = teloTargetIndexGem,
            adaptersFa = adaptersFa,
            conpairPileupNormalPreexist = conpairPileupNormalPreexist,
            conpairPileupTumorPreexist = conpairPileupTumorPreexist,
            serviceAccountKey = serviceAccountKey,
            gcpProject = gcpProject
    }
    # wrap in QC test
    if (costAndQcPass) {
        # swap in dummy file if MultipleMetrics not run (NOTE: in this case you should fill in the in MedianInsertSizes!!)
        File tumorInsertSizeMetrics = select_first([tumorMultipleMetrics.insertSizeMetrics, tumorInsertSizeMetricsPreexist, referenceFa])
        File normalInsertSizeMetrics = select_first([normalMultipleMetrics.insertSizeMetrics, normalInsertSizeMetricsPreexist, referenceFa])
        call baf.NygcBaf {
            input:
                referenceFa = referenceFa,
                local = local,
                pairInfo = pairInfo,
                haplotypecallerAnnotatedVcf = select_first([NygcCancerGermline.haplotypecallerAnnotatedVcf, haplotypecallerAnnotatedVcfPreexist]),
                alleleCountsTxtPreexist = alleleCountsTxtPreexist,
                serviceAccountKey = serviceAccountKey,
                gcpProject = gcpProject
        }
        call calling.NygcCalling {
            input:
                pairInfo = pairInfo,
                library = library,
                intervalListBed = intervalListBed,
                mutect2ChrXRawVcfs = mutect2ChrXRawVcfs,
                mutect2ChrXRawStats = mutect2ChrXRawStats,
                referenceFa = referenceFa,
                tumorInsertSizeMetrics = tumorInsertSizeMetrics,
                normalInsertSizeMetrics = normalInsertSizeMetrics,
                normalMedianInsertSize = normalMedianInsertSize,
                tumorMedianInsertSize = tumorMedianInsertSize,
                local = local,
                listOfChroms = listOfChroms,
                listOfChromsFull = listOfChromsFull,
                callerIntervals = callerIntervals,
                callerIntervalsBedFewNodes = callerIntervalsBedFewNodes,
                invertedIntervalListBed = invertedIntervalListBed,
                callRegions = callRegions,
                splitBedsWgs = splitBedsWgs,
                splitBedsWgsFewNodes = splitBedsWgsFewNodes,
                chromBeds = chromBeds,
                readLength = readLength,
                coordReadLength = coordReadLength,
                uniqCoords = uniqCoords,
                bicseq2ConfigFile = bicseq2ConfigFile,
                bicseq2SegConfigFile = bicseq2SegConfigFile,
                lambda = lambda,
                chromFastas = chromFastas,
                bwaReference = bwaReference,
                bsGenome = bsGenome,
                ponTarGz = ponTarGz,
                gridssAdditionalReference = gridssAdditionalReference,
                configureStrelkaSomaticWorkflow = configureStrelkaSomaticWorkflow,
                lancetJsonLog = lancetJsonLog,
                mantaJsonLog = mantaJsonLog,
                strelkaJsonLog = strelkaJsonLog,
                mutectJsonLog = mutectJsonLog,
                mutectJsonLogFilter = mutectJsonLogFilter,
                gridssPreMemoryGb = gridssPreMemoryGb,
                gridssFilterMemoryGb = gridssFilterMemoryGb,
                gridssHighMem = gridssHighMem,
                mantaHighMem = mantaHighMem,
                mutect2HighMem = mutect2HighMem,
                mutect2Preexist = mutect2Preexist,
                candidateSmallIndelsPreexist = candidateSmallIndelsPreexist,
                diploidSVPreexist = diploidSVPreexist,
                somaticSVPreexist = somaticSVPreexist,
                candidateSVPreexist = candidateSVPreexist,
                unfilteredMantaSVPreexist = unfilteredMantaSVPreexist,
                filteredMantaSVPreexist = filteredMantaSVPreexist,
                strelka2SnvsPreexist = strelka2SnvsPreexist,
                strelka2IndelsPreexist = strelka2IndelsPreexist,
                strelka2SnvPreexist = strelka2SnvPreexist,
                strelka2IndelPreexist = strelka2IndelPreexist,
                lancetPreexist = lancetPreexist,
                gridssVcfPreexist = gridssVcfPreexist,
                gridssUnfilteredVcfChromsPreexist = gridssUnfilteredVcfChromsPreexist,
                bicseq2PngPreexist = bicseq2PngPreexist,
                bicseq2Preexist = bicseq2Preexist,
                serviceAccountKey = serviceAccountKey,
                gcpProject = gcpProject
        }
        
    call postProcessing.NygcPostProcessing {
        input:
            external = external,
            production = production,
            local = local,
            cytoBand = cytoBand,
            dgv = dgv,
            thousandG = thousandG,
            cosmicUniqueBed = cosmicUniqueBed,
            cancerCensusBed = cancerCensusBed,
            ensemblUniqueBed = ensemblUniqueBed,
            gap = gap,
            dgvBedpe = dgvBedpe,
            thousandGVcf = thousandGVcf,
            svPon = svPon,
            cosmicBedPe = cosmicBedPe,
            annotateBicSeq2CnvMem = annotateBicSeq2CnvMem,
            cosmicSigs = cosmicSigs,
            cosmicSbsReference = cosmicSbsReference,
            cosmicDbsReference = cosmicDbsReference,
            cosmicIdReference = cosmicIdReference,
            ponWGSFile = ponWGSFile,
            ponExomeFile = ponExomeFile,
            germFile = germFile,
            referenceFa = referenceFa,
            listOfChroms = listOfChroms,
            deepIntronicsVcf = deepIntronicsVcf,
            clinvarIntronicsVcf = clinvarIntronicsVcf,
            vepGenomeBuild = vepGenomeBuild,
            cosmicCodingChrom = cosmicCodingChrom,
            cosmicNoncodingChrom = cosmicNoncodingChrom,
            dbNSFPReplacementLogic = dbNSFPReplacementLogic,
            MaxEntScan = MaxEntScan,
            dbNSFP4Chrom = dbNSFP4Chrom,
            dbscSNVChrom = dbscSNVChrom,
            gnomadGenomes = gnomadGenomes,
            gnomadExomes = gnomadExomes,
            intervalListBed = intervalListBed,
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
            pairInfo = pairInfo,
            mutect2 = NygcCalling.mutect2,
            filteredMantaSV = NygcCalling.filteredMantaSV,
            strelka2Snv = NygcCalling.strelka2Snv,
            strelka2Indel = NygcCalling.strelka2Indel,
            lancet = NygcCalling.lancet,
            gridssVcf = NygcCalling.gridssVcf,
            bicseq2 = NygcCalling.bicseq2,
            mergedVcfPreexist = mergedVcfPreexist,
            mainVcfPreexist = mainVcfPreexist,
            supplementalVcfPreexist = supplementalVcfPreexist,
            vcfAnnotatedTxtPreexist = vcfAnnotatedTxtPreexist,
            mafPreexist = mafPreexist,
            cnvAnnotatedFinalBedPreexist = cnvAnnotatedFinalBedPreexist,
            cnvAnnotatedSupplementalBedPreexist = cnvAnnotatedSupplementalBedPreexist,
            svFinalBedPePreexist = svFinalBedPePreexist,
            svHighConfidenceFinalBedPePreexist = svHighConfidenceFinalBedPePreexist,
            svSupplementalBedPePreexist = svSupplementalBedPePreexist,
            svHighConfidenceSupplementalBedPePreexist = svHighConfidenceSupplementalBedPePreexist,
            sigsPreexist = sigsPreexist,
            countsPreexist = countsPreexist,
            sig_inputPreexist = sig_inputPreexist,
            reconstructedPreexist = reconstructedPreexist,
            diffPreexist = diffPreexist,
            filteredVcfPreexist = filteredVcfPreexist,
            outputSbsMatrixPreexist = outputSbsMatrixPreexist,
            outputDbsMatrixPreexist = outputDbsMatrixPreexist,
            outputIdMatrixPreexist = outputIdMatrixPreexist,
            sigsSbsPreexist = sigsSbsPreexist,
            sigsDbsPreexist = sigsDbsPreexist,
            sigsIdPreexist = sigsIdPreexist,
            serviceAccountKey = serviceAccountKey,
            gcpProject = gcpProject
    }
    # create resources to replace full BAMs for most review
    call utils.FineGrainedCov as tumorFineGrainedCov {
        input:
            finalBam = pairInfo.tumorFinalBam,
            sampleId = pairInfo.tumorId,
            analysisId = pairInfo.analysisTumorId,
            diskSize = tumorDiskSize + 50
    }
    call utils.FineGrainedCov as normalFineGrainedCov {
        input:
            finalBam = pairInfo.normalFinalBam,
            sampleId = pairInfo.normalId,
            analysisId = pairInfo.analysisNormalId,
            diskSize = normalDiskSize + 50
    }
    String tumorVariantCramPath = "~{pairInfo.analysisTumorId}.variantRegions.cram"
    String tumorVariantCraiPath = "~{pairInfo.analysisTumorId}.variantRegions.crai"
    call utils.MakeVariantCram as tumorMakeVariantCram {
        input:
            finalBam = pairInfo.tumorFinalBam,
            sampleId = pairInfo.tumorId,
            mainVcf = NygcPostProcessing.mainVcf,
            cnvAnnotatedFinalBed = select_first([NygcPostProcessing.cnvAnnotatedFinalBed]),
            svFinalBedPe = select_first([NygcPostProcessing.svFinalBedPe]),
            referenceFa = referenceFa,
            variantCramPath = tumorVariantCramPath,
            variantCraiPath = tumorVariantCraiPath,
            serviceAccountKey = serviceAccountKey,
            gcpProject = gcpProject,
            diskSize = tumorDiskSize + 50
    }
    String normalVariantCramPath = "~{pairInfo.analysisNormalId}.variantRegions.cram"
    String normalVariantCraiPath = "~{pairInfo.analysisNormalId}.variantRegions.crai"
    call utils.MakeVariantCram as normalMakeVariantCram {
        input:
            finalBam = pairInfo.normalFinalBam,
            sampleId = pairInfo.normalId,
            mainVcf = NygcPostProcessing.mainVcf,
            cnvAnnotatedFinalBed = select_first([NygcPostProcessing.cnvAnnotatedFinalBed]),
            svFinalBedPe = select_first([NygcPostProcessing.svFinalBedPe]),
            referenceFa = referenceFa,
            variantCramPath = normalVariantCramPath,
            variantCraiPath = normalVariantCraiPath,
            serviceAccountKey = serviceAccountKey,
            gcpProject = gcpProject,
            diskSize = normalDiskSize + 50
    }
    FinalVcfPairInfo finalVcfPairInfo = object {
                pairId : pairInfo.pairId,
                tumor : pairInfo.tumorId,
                normal : pairInfo.normalId,
                mainVcf : NygcPostProcessing.mainVcf,
                supplementalVcf : NygcPostProcessing.supplementalVcf,
                vcfAnnotatedTxt : NygcPostProcessing.vcfAnnotatedTxt,
                maf : NygcPostProcessing.maf,
                filteredMantaSV : NygcCalling.filteredMantaSV,
                strelka2Snv : NygcCalling.strelka2Snv,
                strelka2Indel : NygcCalling.strelka2Indel,
                mutect2 : NygcCalling.mutect2,
                lancet : NygcCalling.lancet,
                gridssVcf : NygcCalling.gridssVcf,
                bicseq2Png : NygcCalling.bicseq2Png,
                bicseq2 : NygcCalling.bicseq2,
                cnvAnnotatedFinalBed : NygcPostProcessing.cnvAnnotatedFinalBed,
                cnvAnnotatedSupplementalBed : NygcPostProcessing.cnvAnnotatedSupplementalBed,
                svFinalBedPe : NygcPostProcessing.svFinalBedPe,
                svHighConfidenceFinalBedPe : NygcPostProcessing.svHighConfidenceFinalBedPe,
                svSupplementalBedPe : NygcPostProcessing.svSupplementalBedPe,
                svHighConfidenceSupplementalBedPe : NygcPostProcessing.svHighConfidenceSupplementalBedPe,
            }
    } # end of QC test wrap

    output {
        Boolean runAssumeCallingSpeedyOutput = runAssumeCallingSpeedy
        Boolean? costAndQcCheckOutput = costAndQcCheck
        Boolean somaticQcPassOutput = somaticQcPass
        Boolean allCostPass = costPass
        File? mutect2chrXRateTableLog = Mutect2RateTest.mutect2chrXRateTable
        File? mutect2Rate = Mutect2RateTest.mutect2Rate
        File? tumorDistMosdepth = tumorFineGrainedCov.distMosdepthRegion
        File? tumorSummaryMosdepth = tumorFineGrainedCov.summaryMosdepth
        File? tumorRegionsMosdepth = tumorFineGrainedCov.regionsMosdepth
        File? tumorRegionsMosdepthIndex = tumorFineGrainedCov.regionsMosdepthIndex
        File? normalDistMosdepth = normalFineGrainedCov.distMosdepthRegion
        File? normalSummaryMosdepth = normalFineGrainedCov.summaryMosdepth
        File? normalRegionsMosdepth = normalFineGrainedCov.regionsMosdepth
        File? normalRegionsMosdepthIndex = normalFineGrainedCov.regionsMosdepthIndex
        Cram? tumorVariantCram = tumorMakeVariantCram.variantCram
        Cram? normalVariantCram = normalMakeVariantCram.variantCram
        FinalVcfPairInfo? finalVcfPairInfos = finalVcfPairInfo
        Array[File?] alignmentSummaryMetrics = select_all([tumorAlignmentSummaryMetricsPreexist, normalAlignmentSummaryMetricsPreexist,
                                                         tumorMultipleMetrics.alignmentSummaryMetrics, normalMultipleMetrics.alignmentSummaryMetrics])
        Array[File?] qualityByCyclePdf = select_all([tumorQualityByCyclePdfPreexist, normalQualityByCyclePdfPreexist,
                                                        tumorMultipleMetrics.qualityByCyclePdf, normalMultipleMetrics.qualityByCyclePdf])
        Array[File?] baseDistributionByCycleMetrics = select_all([tumorBaseDistributionByCycleMetricsPreexist, normalBaseDistributionByCycleMetricsPreexist,
                                                        tumorMultipleMetrics.baseDistributionByCycleMetrics, normalMultipleMetrics.baseDistributionByCycleMetrics])
        Array[File?] qualityByCycleMetrics = select_all([tumorQualityByCycleMetricsPreexist, normalQualityByCycleMetricsPreexist,
                                                        tumorMultipleMetrics.qualityByCycleMetrics, normalMultipleMetrics.qualityByCycleMetrics])
        Array[File?] baseDistributionByCyclePdf = select_all([tumorBaseDistributionByCyclePdfPreexist, normalBaseDistributionByCyclePdfPreexist,
                                                            tumorMultipleMetrics.baseDistributionByCyclePdf, normalMultipleMetrics.baseDistributionByCyclePdf])
        Array[File?] qualityDistributionPdf = select_all([tumorQualityDistributionPdfPreexist, normalQualityDistributionPdfPreexist,
                                                        tumorMultipleMetrics.qualityDistributionPdf, normalMultipleMetrics.qualityDistributionPdf])
        Array[File?] qualityDistributionMetrics = select_all([tumorQualityDistributionMetricsPreexist, normalQualityDistributionMetricsPreexist,
                                                                tumorMultipleMetrics.qualityDistributionMetrics, normalMultipleMetrics.qualityDistributionMetrics])
        Array[File?] insertSizeHistogramPdf = select_all([tumorInsertSizeHistogramPdfPreexist, normalInsertSizeHistogramPdfPreexist,
                                                        tumorMultipleMetrics.insertSizeHistogramPdf, normalMultipleMetrics.insertSizeHistogramPdf])
        Array[File?] insertSizeMetrics = select_all([tumorInsertSizeMetrics, normalInsertSizeMetrics,
                                                        tumorMultipleMetrics.insertSizeMetrics, normalMultipleMetrics.insertSizeMetrics])
#         # Mutect2
#         File mutect2 = NygcCalling.mutect2
#         # Manta
#         File filteredMantaSV = NygcCalling.filteredMantaSV
#         # Strelka2
#         File strelka2Snv = NygcCalling.strelka2Snv
#         File strelka2Indel = NygcCalling.strelka2Indel
#         # Lancet
#         File lancet = NygcCalling.lancet
#         # Gridss
#         IndexedVcf gridssVcf = NygcCalling.gridssVcf
#         # Bicseq2
#         File? bicseq2Png = NygcCalling.bicseq2Png
#         File? bicseq2 = NygcCalling.bicseq2
        IndexedVcf? haplotypecallerFinalFiltered = select_first([NygcCancerGermline.haplotypecallerFinalFiltered, haplotypecallerFinalFilteredPreexist])
        File? filteredHaplotypecallerAnnotatedVcf = select_first([NygcCancerGermline.filteredHaplotypecallerAnnotatedVcf, filteredHaplotypecallerAnnotatedVcfPreexist])
        File? haplotypecallerAnnotatedVcf = select_first([NygcCancerGermline.haplotypecallerAnnotatedVcf, haplotypecallerAnnotatedVcfPreexist])
        File? kouramiResult = select_first([NygcCancerGermline.kouramiResult, kouramiResultPreexist])
        File? r1HlaFastq = select_first([NygcCancerGermline.r1HlaFastq, r1HlaFastqPreexist])
        File? r2HlaFastq = select_first([NygcCancerGermline.r2HlaFastq, r2HlaFastqPreexist])
        # ancestry
        File? beagleFileContinental = select_first([NygcCancerGermline.beagleFileContinental, beagleFileContinentalPreexist])
        File? fastNgsAdmixQoptContinental = select_first([NygcCancerGermline.fastNgsAdmixQoptContinental, fastNgsAdmixQoptContinentalPreexist])
        File? beagleFilePopulation = select_first([NygcCancerGermline.beagleFilePopulation, beagleFilePopulationPreexist])
        File? fastNgsAdmixQoptPopulation = select_first([NygcCancerGermline.fastNgsAdmixQoptPopulation, fastNgsAdmixQoptPopulationPreexist])
        File? alleleCountsTxt = NygcBaf.alleleCountsTxt
        # Conpair
        File concordanceAll = NygcPairedTasks.concordanceAll
        File concordanceHomoz = NygcPairedTasks.concordanceHomoz
        File contamination = NygcPairedTasks.contamination
        File conpairPileupTumor = NygcPairedTasks.conpairPileupTumor
        File conpairPileupNormal = NygcPairedTasks.conpairPileupNormal
        # TeloRead Descriptions
        File normalGemLength = NygcPairedTasks.normalGemLength
        # Aligned/Unaligned Telomeric reads
        File normalTeloR1Fastq = NygcPairedTasks.normalTeloR1Fastq
        File normalTeloR2Fastq = NygcPairedTasks.normalTeloR2Fastq
        File normalSingMappedFastq = NygcPairedTasks.normalSingMappedFastq
        File normalAlignmentCsv = NygcPairedTasks.normalAlignmentCsv
        # normal teloreads
        # TelFusDetector
        File normalFusionsFiltered = NygcPairedTasks.normalFusionsFiltered
        File normalFusionsPass = NygcPairedTasks.normalFusionsPass
        File normalAllChromosomesCov = NygcPairedTasks.normalAllChromosomesCov
        File normalFusionsRates = NygcPairedTasks.normalFusionsRates
        # telomeasures
        File telomeasuresNormal = NygcPairedTasks.telomeasuresNormal
        # tumor teloreads
        # TeloRead Descriptions
        File tumorGemLength = NygcPairedTasks.tumorGemLength
        # Aligned/Unaligned Telomeric reads
        File tumorTeloR1Fastq = NygcPairedTasks.tumorTeloR1Fastq
        File tumorTeloR2Fastq = NygcPairedTasks.tumorTeloR2Fastq
        File tumorSingMappedFastq = NygcPairedTasks.tumorSingMappedFastq
        File tumorAlignmentCsv = NygcPairedTasks.tumorAlignmentCsv
        # TelFusDetector
        File tumorFusionsFiltered = NygcPairedTasks.tumorFusionsFiltered
        File tumorFusionsPass = NygcPairedTasks.tumorFusionsPass
        File tumorAllChromosomesCov = NygcPairedTasks.normalAllChromosomesCov
        File tumorFusionsRates = NygcPairedTasks.tumorFusionsRates
        # telomeasures
        File telomeasuresTumor = NygcPairedTasks.telomeasuresTumor
        # MSI
        File mantisWxsKmerCountsFinal = NygcPairedTasks.mantisWxsKmerCountsFinal
        File mantisWxsKmerCountsFiltered = NygcPairedTasks.mantisWxsKmerCountsFiltered
        File mantisExomeTxt = NygcPairedTasks.mantisExomeTxt
        File mantisStatusFinal = NygcPairedTasks.mantisStatusFinal
        # -> Post processing
#         File mainVcf = NygcPostProcessing.mainVcf
        File? filteredVcf = NygcPostProcessing.filteredVcf # highConfidenceVcf
#         File supplementalVcf = NygcPostProcessing.supplementalVcf
#         File vcfAnnotatedTxt = NygcPostProcessing.vcfAnnotatedTxt
#         File maf = NygcPostProcessing.maf
#         File? cnvAnnotatedFinalBed  = NygcPostProcessing.cnvAnnotatedFinalBed
#         File? cnvAnnotatedSupplementalBed = NygcPostProcessing.cnvAnnotatedSupplementalBed
#         File? svFinalBedPe = NygcPostProcessing.svFinalBedPe
#         File? svHighConfidenceFinalBedPe = NygcPostProcessing.svHighConfidenceFinalBedPe
#         File? svSupplementalBedPe = NygcPostProcessing.svSupplementalBedPe
#         File? svHighConfidenceSupplementalBedPe = NygcPostProcessing.svHighConfidenceSupplementalBedPe
        # sigs
        File? sigs = NygcPostProcessing.sigs
        File? counts = NygcPostProcessing.counts
        File? sig_input = NygcPostProcessing.sig_input
        File? reconstructed = NygcPostProcessing.reconstructed
        File? diff = NygcPostProcessing.diff
        File? outputSbsMatrix = NygcPostProcessing.outputSbsMatrix
        File? outputDbsMatrix = NygcPostProcessing.outputDbsMatrix
        File? outputIdMatrix = NygcPostProcessing.outputIdMatrix
        File? sigsSbs = NygcPostProcessing.sigsSbs
        Array[File]? sigsDbs = NygcPostProcessing.sigsDbs
        Array[File]? sigsId = NygcPostProcessing.sigsId
        
    }
}
