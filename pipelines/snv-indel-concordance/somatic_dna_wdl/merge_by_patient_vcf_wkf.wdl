version 1.0

import "merge_vcf/merge_vcf_participant_wkf.wdl" as mergeParticipantVcf
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

workflow MergeParticipantVcfs {
    # command
    #   Call variants in BAMs
    #   merge and filter merged VCFs
    #   annotate
    input {
        Boolean local = false
        Array[File] centerVcfs
        Array[String] tumorBarcodeAliquots
        Array[String] nonNormalKinds
        String normalBarcodeAliquot
        Array[Bam] tumorFinalBams
        Bam normalFinalBam
        String participantId
        
        # merge callers
        File intervalListBed

        # annotation:
        String vepGenomeBuild
        Map[String, IndexedVcf] cosmicCodingChrom
        Map[String, IndexedVcf] cosmicNoncodingChrom
        Array[String]+ listOfChroms
        File cosmicCensus

        # Public
        File cancerResistanceMutations


        # post annotation
        File cosmicCensus

        IndexedReference referenceFa
        # optional annotation
        Boolean reAnnotate = false
        String vepGenomeBuild
        Map[String, IndexedVcf] cosmicCodingChrom
        Map[String, IndexedVcf] cosmicNoncodingChrom
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
        IndexedVcf deepIntronicsVcf
        IndexedVcf clinvarIntronicsVcf
        
        File? serviceAccountKey
        String? gcpProject
    }
    
    call mergeParticipantVcf.MergeParticipantVcf {
        input:
            local = local,
            centerVcfs = centerVcfs,
            tumorBarcodeAliquots = tumorBarcodeAliquots,
            nonNormalKinds = nonNormalKinds,
            normalBarcodeAliquot = normalBarcodeAliquot,
            tumorFinalBams = tumorFinalBams,
            normalFinalBam = normalFinalBam,
            participantId = participantId,
            listOfChroms = listOfChroms,
            referenceFa = referenceFa,
            cosmicCensus = cosmicCensus,
            
            reAnnotate = reAnnotate,
            vepGenomeBuild = vepGenomeBuild,
            cosmicCodingChrom = cosmicCodingChrom,
            cosmicNoncodingChrom = cosmicNoncodingChrom,
            vepCacheChrom = vepCacheChrom,
            annotationsChrom = annotationsChrom,
            plugins = plugins,
            dbNSFPReplacementLogic = dbNSFPReplacementLogic,
            MaxEntScan = MaxEntScan,
            dbNSFP4Chrom = dbNSFP4Chrom,
            dbscSNVChrom = dbscSNVChrom,
            gnomadGenomes = gnomadGenomes,
            gnomadExomes = gnomadExomes,
            vepFastaReference = vepFastaReference,
            cancerResistanceMutations = cancerResistanceMutations,
            deepIntronicsVcf = deepIntronicsVcf,
            clinvarIntronicsVcf = clinvarIntronicsVcf,
            
            serviceAccountKey = serviceAccountKey,
            gcpProject = gcpProject
    }
    
    output {
        Array[File] unionPairedVcf = MergeParticipantVcf.unionPairedVcf
        Array[File] intersectionPairedVcf = MergeParticipantVcf.intersectionPairedVcf 
        File unionVcf = MergeParticipantVcf.unionVcf
        File intersectionVcf = MergeParticipantVcf.intersectionVcf
        File mutect2Genotyped = MergeParticipantVcf.mutect2Genotyped
        File mutect2UnfilteredGenotyped = MergeParticipantVcf.mutect2UnfilteredGenotyped
    }

    
    
}