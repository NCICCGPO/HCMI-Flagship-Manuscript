version 1.0

import "../wdl_structs.wdl"
import "../merge_vcf/merge_vcf_wkf.wdl" as mergeVcf
import "../annotate/annotate_wkf.wdl" as annotate
import "../annotate/annotate_cnv_sv_wkf.wdl" as annotate_cnv_sv
import "../variant_analysis/deconstruct_sigs_wkf.wdl" as deconstructSigs
import "../variant_analysis/musical_sigs_wkf.wdl" as musicalSigs
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


workflow NygcPostProcessing {
    input {
        Boolean external = false
        Boolean local = false
        Boolean production = true
        IndexedReference referenceFa
        Array[String]+ listOfChroms
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
        # merge callers
        File intervalListBed
        String library
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
        # annotate SNV/INDELs
        String vepGenomeBuild
        Map[String, IndexedVcf] cosmicCodingChrom
        Map[String, IndexedVcf] cosmicNoncodingChrom
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
        File ensemblEntrez
        Int annotateBicSeq2CnvMem = 36
        # patient-specific inputs
        PairInfo pairInfo
        # Mutect2
        File mutect2
        # Manta
        File filteredMantaSV
        # Strelka2
        File strelka2Snv
        File strelka2Indel
        # Lancet
        File lancet 
        # Gridss
        IndexedVcf gridssVcf
        # Bicseq2
        File? bicseq2
        # needed for work on cloud without localizing
        File? serviceAccountKey
        String? gcpProject
        # Pre-existing finished files (in case only some steps need to rerun)
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
    }
    if ((!defined(mergedVcfPreexist)) && (!defined(mainVcfPreexist))) {
        PreMergedPairVcfInfo preMergedPairVcfInfo = object {
            pairId : pairInfo.analysisPairId,
            filteredMantaSV : filteredMantaSV,
            strelka2Snv : strelka2Snv,
            strelka2Indel : strelka2Indel,
            mutect2 : mutect2,
            lancet : lancet,
            tumor : pairInfo.tumorId,
            normal : pairInfo.normalId,
            tumorFinalBam : pairInfo.tumorFinalBam,
            normalFinalBam : pairInfo.normalFinalBam
        }
        if (library == 'WGS') {
            call mergeVcf.MergeVcf as wgsMergeVcf {
                input:
                    external = external,
                    preMergedPairVcfInfo = preMergedPairVcfInfo,
                    referenceFa = referenceFa,
                    listOfChroms = listOfChroms,
                    intervalListBed = intervalListBed,
                    ponFile = ponWGSFile,
                    library = library,
                    serviceAccountKey = serviceAccountKey,
                    gcpProject = gcpProject
            }
        }
        if (library == 'Exome') {
            call mergeVcf.MergeVcf as exomeMergeVcf {
                input:
                    external = external,
                    preMergedPairVcfInfo = preMergedPairVcfInfo,
                    referenceFa = referenceFa,
                    listOfChroms = listOfChroms,
                    intervalListBed = intervalListBed,
                    ponFile = ponExomeFile,
                    library = library,
                    serviceAccountKey = serviceAccountKey,
                    gcpProject = gcpProject
            }
        }
    }
    if (!defined(mainVcfPreexist)) {
        File mergedVcf = select_first([wgsMergeVcf.mergedVcf, exomeMergeVcf.mergedVcf, mergedVcfPreexist])
        call annotate.Annotate {
            input:
                unannotatedVcf = mergedVcf,
                referenceFa = referenceFa,
                production = production,
                local = local,
                tumorId = pairInfo.tumorId,
                normalId = pairInfo.normalId,
                pairName = pairInfo.analysisPairId,
                vepGenomeBuild = vepGenomeBuild,
                cosmicCodingChrom = cosmicCodingChrom,
                cosmicNoncodingChrom = cosmicNoncodingChrom,
                dbNSFPReplacementLogic = dbNSFPReplacementLogic,
                MaxEntScan = MaxEntScan,
                dbNSFP4Chrom = dbNSFP4Chrom,
                dbscSNVChrom = dbscSNVChrom,
                gnomadGenomes = gnomadGenomes,
                gnomadExomes = gnomadExomes,
                listOfChroms = listOfChroms,
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
        File mainVcfFinal = Annotate.pairVcfInfo.mainVcf
        File supplementalVcfFinal = Annotate.pairVcfInfo.supplementalVcf
        File vcfAnnotatedTxtFinal = Annotate.pairVcfInfo.vcfAnnotatedTxt
        File mafFinal = Annotate.pairVcfInfo.maf
    }
    if (defined(bicseq2)){
        File bicseq2Run = select_first([bicseq2, cytoBand])
        if (!defined(cnvAnnotatedFinalBedPreexist)){
            call annotate_cnv_sv.AnnotateCnvSv {
                input:
                    tumor=pairInfo.tumorId,
                    normal=pairInfo.normalId,
                    pairName=pairInfo.analysisPairId,
                    listOfChroms=listOfChroms,
                    bicseq2=bicseq2Run,
                    cytoBand=cytoBand,
                    dgv=dgv,
                    thousandG=thousandG,
                    cosmicUniqueBed=cosmicUniqueBed,
                    cancerCensusBed=cancerCensusBed,
                    ensemblUniqueBed=ensemblUniqueBed,
                    filteredMantaSV=filteredMantaSV,
                    gridssVcf=gridssVcf,
                    vepGenomeBuild=vepGenomeBuild,
                    gap=gap,
                    dgvBedpe=dgvBedpe,
                    thousandGVcf=thousandGVcf,
                    svPon=svPon,
                    cosmicBedPe=cosmicBedPe,
                    annotateBicSeq2CnvMem = annotateBicSeq2CnvMem
            }
        }
    }
    File mainVcfRun = select_first([mainVcfFinal, mainVcfPreexist])
    if (!defined(filteredVcfPreexist)){
        call tests.FilterHighConfidence {
            input:
                vcf = mainVcfRun,
                filteredVcfPath = "~{pairInfo.pairId}.snv.indel.high_confidence.v7.annotated.vcf",
        }
    }
    File filteredVcfFinal = select_first([FilterHighConfidence.filteredVcf, filteredVcfPreexist])
    # these sections are quick to rerun
    if (!defined(sigsPreexist)){
        call deconstructSigs.DeconstructSig {
            input:
                pairId = pairInfo.analysisPairId,
                mainVcf = filteredVcfFinal,
                cosmicSigs = cosmicSigs,
                vepGenomeBuild = vepGenomeBuild
        }
    }
    File sigsFinal = select_first([DeconstructSig.sigs, sigsPreexist])
    File countsFinal = select_first([DeconstructSig.counts, countsPreexist])
    File sig_inputFinal = select_first([DeconstructSig.sigInput, sig_inputPreexist])
    File reconstructedFinal = select_first([DeconstructSig.reconstructed, reconstructedPreexist])
    File diffFinal = select_first([DeconstructSig.diff, diffPreexist])
    if (!defined(outputSbsMatrixPreexist)){
        call musicalSigs.MusicalSig {
            input:
                pairId = pairInfo.analysisPairId,
                mainVcf = filteredVcfFinal,
                vepGenomeBuild = vepGenomeBuild,
                cosmicSbsReference = cosmicSbsReference,
                cosmicDbsReference = cosmicDbsReference,
                cosmicIdReference = cosmicIdReference
        }
    }
    File outputSbsMatrixFinal = select_first([MusicalSig.outputSbsMatrix, outputSbsMatrixPreexist])
    File outputDbsMatrixFinal = select_first([MusicalSig.outputDbsMatrix, outputDbsMatrixPreexist])
    File outputIdMatrixFinal = select_first([MusicalSig.outputIdMatrix, outputIdMatrixPreexist])
    File sigsSbsFinal = select_first([MusicalSig.sigsSbs, sigsSbsPreexist])
    Array[File] sigsDbsFinal = select_all([MusicalSig.sigsDbs, sigsDbsPreexist])
    Array[File] sigsIdFinal = select_all([MusicalSig.sigsId, sigsIdPreexist])
    output {
        File mainVcf = mainVcfRun
        File filteredVcf = filteredVcfFinal
        File supplementalVcf = select_first([supplementalVcfFinal, supplementalVcfPreexist])
        File vcfAnnotatedTxt = select_first([vcfAnnotatedTxtFinal, vcfAnnotatedTxtPreexist])
        File maf = select_first([mafFinal, mafPreexist])
        File? cnvAnnotatedFinalBed  = select_first([AnnotateCnvSv.cnvAnnotatedFinalBed, cnvAnnotatedFinalBedPreexist])
        File? cnvAnnotatedSupplementalBed = select_first([AnnotateCnvSv.cnvAnnotatedSupplementalBed,cnvAnnotatedSupplementalBedPreexist])
        File? svFinalBedPe = select_first([AnnotateCnvSv.svFinalBedPe, svFinalBedPePreexist])
        File? svHighConfidenceFinalBedPe = select_first([AnnotateCnvSv.svHighConfidenceFinalBedPe,svHighConfidenceFinalBedPePreexist])
        File? svSupplementalBedPe = select_first([AnnotateCnvSv.svSupplementalBedPe, svSupplementalBedPePreexist])
        File? svHighConfidenceSupplementalBedPe = select_first([AnnotateCnvSv.svHighConfidenceSupplementalBedPe, svHighConfidenceSupplementalBedPePreexist])
        # sigs
        File sigs = sigsFinal
        File counts = countsFinal
        File sig_input = sig_inputFinal
        File reconstructed = reconstructedFinal
        File diff = diffFinal
        File outputSbsMatrix = outputSbsMatrixFinal
        File outputDbsMatrix = outputDbsMatrixFinal
        File outputIdMatrix = outputIdMatrixFinal
        File sigsSbs = sigsSbsFinal
        Array[File] sigsDbs = sigsDbsFinal
        Array[File] sigsId = sigsIdFinal
    }
}