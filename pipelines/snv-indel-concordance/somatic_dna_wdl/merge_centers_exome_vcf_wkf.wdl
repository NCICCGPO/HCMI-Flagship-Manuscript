version 1.0
import "somatic_bam_snv_indel_wkf.wdl" as somaticBamExome
import "test/tests.wdl" as tests
import "calling/calling.wdl" as calling
# import "test/utils.wdl" as utils
import "merge_vcf/merge_only_exome_vcf_centers_wkf.wdl" as mergeVcfExome
import "annotate/annotate_center_wkf.wdl" as annotateVcf
import "wdl_structs.wdl"

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

workflow MergeVcf {
    # command
    #   Call variants in BAMs
    #   merge and filter raw VCFs
    #   annotate
    input {
        Boolean external = false
        File broadMaf
        File nygcVcf
        Bam tumorFinalBam
        Bam normalFinalBam
        String tumorId
        String normalId
        String pairId
        
        # merge callers
        File intervalListBed

        # annotation:
        String vepGenomeBuild
        Map[String, IndexedVcf] cosmicCodingChrom
        Map[String, IndexedVcf] cosmicNoncodingChrom
        Array[String]+ listOfChroms

        # Public
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

        # Public
        File cancerResistanceMutations
        IndexedVcf deepIntronicsVcf
        IndexedVcf clinvarIntronicsVcf

        # post annotation
        File cosmicCensus

        String library = "Exome"

        IndexedReference referenceFa
        IndexedTable exonicRegionBedAnnotationsGz = {"table" : "fake", "index" : "fake"}
        File? serviceAccountKey
        String? gcpProject
        
        
        ##########################
        BwaReference bwaReference
        Boolean production = true
        # For Tumor-Normal QC
        File markerBedFile
        File markerTxtFile

        # calling
        Array[String]+ listOfChromsFull
        Array[String]+ callerIntervals
        Array[String]+ exomeCallerIntervals = callerIntervals
        File invertedIntervalListBed
        IndexedTable callRegions
        Array[File] splitBedsWgs
        Map[String, File] chromBeds
        File lancetJsonLog
        File mantaJsonLog
        File strelkaJsonLog
        File mutectJsonLog
        File mutectJsonLogFilter
        File configureStrelkaSomaticWorkflow

        File ponWGSFile
        File ponExomeFile
        IndexedVcf gnomadBiallelic

        # mantis
        File mantisBed

        # post annotation
        File ensemblEntrez
        # germline
        File excludeIntervalList
        Array[File] scatterIntervalsHcs

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

        # signatures
        File cosmicSigs

        # Gridss resources need a lot of fine grained control
        Int gridssPreMemoryGb = 60
        Int gridssFilterMemoryGb = 32
        Boolean gridssHighMem = false
        Boolean mantaHighMem = false
        Boolean mutect2HighMem = false
        
        # need to create docker with gcloud
        File? serviceAccountKey
        String? gcpProject
        File awsConfig
        File awsCredentials
        String endpointUrl = "https://gdc-jamboree-objstore.datacommons.io"
    }

    call mergeVcfExome.MergeVcfCenters {
        input:
            external = external,
            broadMaf = broadMaf,
            nygcVcf = nygcVcf,
            normalFinalBam=normalFinalBam,
            tumorFinalBam=tumorFinalBam,
            tumorId=tumorId,
            normalId=normalId,
            library=library,
            pairId=pairId,
            referenceFa = referenceFa,
            listOfChroms = listOfChroms,
            intervalListBed = intervalListBed,
            serviceAccountKey = serviceAccountKey,
            gcpProject = gcpProject
    }
    
    call annotateVcf.Annotate {
        input:
            unannotatedVcf = MergeVcfCenters.mergedVcf,
            tumorId=tumorId,
            normalId=normalId,
            pairName = pairId,
            referenceFa = referenceFa,
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
            library = library,
            listOfChroms = listOfChroms,
            exonicRegionBedAnnotationsGz = exonicRegionBedAnnotationsGz,
    
    }
    
    output {
        File mergedVcf = MergeVcfCenters.mergedVcf
        File mainVcf = Annotate.mainVcf
        File unionMaf = Annotate.unionMaf
        File intersectionMaf = Annotate.intersectionMaf
        File mainIntersectionVcf = Annotate.mainIntersectionVcf
    }
}
