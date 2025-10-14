version 1.0

import "../wdl_structs.wdl"
import "../alignment_analysis/kourami_wfk.wdl" as kourami
import "../alignment_analysis/fastngsadmix_wkf.wdl" as fastNgsAdmix
import "../germline/germline_wkf.wdl" as germline
import "../annotate/germline_annotate_wkf.wdl" as germlineAnnotate
import "../annotate/annotate_wkf.wdl" as annotate

import "../tasks/utils.wdl" as utils
import "../test/tests.wdl" as tests

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


workflow NygcCancerGermline {
    input {
        Boolean production = true
        Boolean local = false
        String library = "WGS"
        # kourami
        BwaReference kouramiReference
        File kouramiFastaGem1Index
        File intervalListBed

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
        IndexedReference referenceFa
        # annotation:
        String vepGenomeBuild
        Map[String, IndexedVcf] cosmicCodingChrom
        Map[String, IndexedVcf] cosmicNoncodingChrom
        Array[String]+ listOfChroms
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
        
        # patient-specific inputs
        String analysisNormalId
        SampleBamInfo normalSampleBamInfo
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
        IndexedVcf? haplotypecallerFinalFilteredPreexist
        
        Boolean germlineHighMem = false
        # need to create docker with gcloud
        File? serviceAccountKey
        String? gcpProject
    }
    if (!defined(kouramiResultPreexist)) {
        call kourami.Kourami {
            input:
                sampleId = normalSampleBamInfo.sampleId,
                analysisId = analysisNormalId,
                kouramiReference = kouramiReference,
                finalBam = normalSampleBamInfo.finalBam,
                kouramiFastaGem1Index = kouramiFastaGem1Index,
                referenceFa = referenceFa,
                serviceAccountKey = serviceAccountKey,
                gcpProject = gcpProject
        }
    }
    File kouramiResultFinal = select_first([Kourami.result, kouramiResultPreexist])
    File r1HlaFastqFinal = select_first([Kourami.r1HlaFastq, r1HlaFastqPreexist])
    File r2HlaFastqFinal = select_first([Kourami.r2HlaFastq, r2HlaFastqPreexist])
    if (!defined(fastNgsAdmixQoptPopulationPreexist)) {
        call fastNgsAdmix.FastNgsAdmix as fastNgsAdmixContinental {
            input:
                normalFinalBam = normalSampleBamInfo.finalBam,
                fastNgsAdmixSites = fastNgsAdmixContinentalSites,
                fastNgsAdmixSitesBin = fastNgsAdmixContinentalSitesBin,
                fastNgsAdmixSitesIdx = fastNgsAdmixContinentalSitesIdx,
                fastNgsAdmixChroms = fastNgsAdmixChroms,
                fastNgsAdmixRef = fastNgsAdmixContinentalRef,
                fastNgsAdmixNind = fastNgsAdmixContinentalNind,
                outprefix = analysisNormalId + "_continental"
        }
    
        call fastNgsAdmix.FastNgsAdmix as fastNgsAdmixPopulation {
            input:
                normalFinalBam = normalSampleBamInfo.finalBam,
                fastNgsAdmixSites = fastNgsAdmixPopulationSites,
                fastNgsAdmixSitesBin = fastNgsAdmixPopulationSitesBin,
                fastNgsAdmixSitesIdx = fastNgsAdmixPopulationSitesIdx,
                fastNgsAdmixChroms = fastNgsAdmixChroms,
                fastNgsAdmixRef = fastNgsAdmixPopulationRef,
                fastNgsAdmixNind = fastNgsAdmixPopulationNind,
                outprefix = analysisNormalId + "population"
        }
    }
    if (!defined(haplotypecallerFinalFilteredPreexist)) {
        call germline.Germline {
            input:
                finalBam = normalSampleBamInfo.finalBam,
                normal = normalSampleBamInfo.sampleId,
                outputPrefix = analysisNormalId,
                referenceFa = referenceFa,
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
                highMem = germlineHighMem
        }
    }
    IndexedVcf haplotypecallerFinalFilteredRun = select_first([Germline.haplotypecallerFinalFiltered, haplotypecallerFinalFilteredPreexist])
    if (!defined(filteredHaplotypecallerAnnotatedVcfPreexist)) {
        String filteredHaplotypecallerAnnotatedVcfPath = "~{analysisNormalId}.haplotypecaller.gatk.annotated.vcf"
        call germlineAnnotate.GermlineAnnotate as filteredGermlineAnnotate {
            input:
                unannotatedVcf = haplotypecallerFinalFilteredRun,
                production = production,
                local = local,
                referenceFa = referenceFa,
                normal = normalSampleBamInfo.sampleId,
                haplotypecallerAnnotatedVcfPath = filteredHaplotypecallerAnnotatedVcfPath,
                listOfChroms = listOfChroms,
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
                library = library
        }
    }
    File filteredHaplotypecallerAnnotatedVcfRun = select_first([filteredGermlineAnnotate.haplotypecallerAnnotatedVcf, filteredHaplotypecallerAnnotatedVcfPreexist])
    if (!defined(haplotypecallerAnnotatedVcfPreexist)) {
        String unfilteredHaplotypecallerAnnotatedVcfPath = "~{analysisNormalId}.haplotypecaller.gatk.annotated.unfiltered.vcf"
        call germlineAnnotate.GermlineAnnotate as unFilteredGermlineAnnotate {
            input:
                unannotatedVcf = select_first([Germline.haplotypecallerVcf, haplotypecallerVcfPreexist]),
                haplotypecallerAnnotatedVcfPath = unfilteredHaplotypecallerAnnotatedVcfPath,
                production = production,
                local = local,
                referenceFa = referenceFa,
                normal = normalSampleBamInfo.sampleId,
                listOfChroms = listOfChroms,
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
                library = library
        }
    }
    File haplotypecallerAnnotatedVcfRun = select_first([unFilteredGermlineAnnotate.haplotypecallerAnnotatedVcf, haplotypecallerAnnotatedVcfPreexist])
    output {
        IndexedVcf haplotypecallerFinalFiltered = haplotypecallerFinalFilteredRun
        File filteredHaplotypecallerAnnotatedVcf = filteredHaplotypecallerAnnotatedVcfRun
        File haplotypecallerAnnotatedVcf = haplotypecallerAnnotatedVcfRun
        File kouramiResult = kouramiResultFinal
        File r1HlaFastq = r1HlaFastqFinal
        File r2HlaFastq = r2HlaFastqFinal
        # ancestry
        File beagleFileContinental = select_first([fastNgsAdmixContinental.beagleFile, beagleFileContinentalPreexist])
        File fastNgsAdmixQoptContinental = select_first([fastNgsAdmixContinental.fastNgsAdmixQopt, fastNgsAdmixQoptContinentalPreexist])
        File beagleFilePopulation = select_first([fastNgsAdmixPopulation.beagleFile, beagleFilePopulationPreexist])
        File fastNgsAdmixQoptPopulation = select_first([fastNgsAdmixPopulation.fastNgsAdmixQopt, fastNgsAdmixQoptPopulationPreexist])
    }
}