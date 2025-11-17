version 1.0

import "merge_vcf/merge_w_exome_vcf_centers_wkf.wdl" as mergeVcfExome
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
        File broadVcf
        File broadExomeVcf
        File nygcVcf
        File washuVcf
        File washuExomeVcf
        Bam normalFinalBam
        Bam tumorFinalBam
        String tumorId
        String normalId
        String tumorIdExome
        String normalIdExome
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

        String library

        IndexedReference referenceFa
        IndexedTable exonicRegionBedAnnotationsGz = {"table" : "fake", "index" : "fake"}
        File? serviceAccountKey
        String? gcpProject
    }
    
    call mergeVcfExome.MergeVcfCenters {
        input:
            external = external,
            broadVcf = broadVcf,
            broadExomeVcf = broadExomeVcf,
            nygcVcf = nygcVcf,
            washuVcf = washuVcf,
            washuExomeVcf = washuExomeVcf,
            normalFinalBam=normalFinalBam,
            tumorFinalBam=tumorFinalBam,
            tumorId=tumorId,
            normalId=normalId,
            tumorIdExome=tumorIdExome,
            normalIdExome=normalIdExome,
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
