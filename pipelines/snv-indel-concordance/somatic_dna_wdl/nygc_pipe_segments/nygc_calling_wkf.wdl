version 1.0

import "../calling/calling.wdl" as callingTasks
import "../wdl_structs.wdl"
import "../calling/mutect2_wkf.wdl" as mutect2
import "../calling/strelka2_wkf.wdl" as strelka2
import "../calling/manta_wkf.wdl" as manta
import "../calling/lancet_wkf.wdl" as lancet
import "../calling/gridss_wkf.wdl" as gridss
import "../calling/bicseq2_wkf.wdl" as bicseq2
import "../tasks/utils.wdl" as utils

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


workflow NygcCalling {
    # command
    #   Call variants in BAMs
    #   merge and filter raw VCFs
    #   annotate
    input {
        Boolean local = false
        String library

        PairInfo pairInfo
        #   mutect2
        Array[String]+ listOfChroms
        Array[String]+ listOfChromsFull
        Array[String]+ callerIntervals
        Array[File]+ callerIntervalsBedFewNodes
        File invertedIntervalListBed
        IndexedReference referenceFa
        Array[File] mutect2ChrXRawVcfs = []
        Array[File] mutect2ChrXRawStats = []
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
        File tumorInsertSizeMetrics
        File normalInsertSizeMetrics
        Int lambda = 4
        # Gridss
        BwaReference bwaReference
        String bsGenome
        File ponTarGz
        Array[File] gridssAdditionalReference
        # Strelka2
        File configureStrelkaSomaticWorkflow
        File intervalListBed

        File lancetJsonLog
        File mantaJsonLog
        File strelkaJsonLog
        File mutectJsonLog
        File mutectJsonLogFilter

        # Gridss resources need a lot of fine grained control
        Int gridssPreMemoryGb = 60
        Int gridssFilterMemoryGb = 32
        Boolean gridssHighMem = false
        Boolean mantaHighMem = false
        Boolean mutect2HighMem = false
        
        # needed for work on cloud without localizing
        File? serviceAccountKey
        String? gcpProject
        
        # finished files (in case only some callers need to rerun)
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

    }
    if (!defined(mutect2Preexist)){
        call mutect2.Mutect2 {
            input:
                local = local,
                library = library,
                invertedIntervalListBed = invertedIntervalListBed,
                callerIntervals = callerIntervals,
                callerIntervalsBedFewNodes = callerIntervalsBedFewNodes,
                mutectJsonLogFilter = mutectJsonLogFilter,
                mutectJsonLog = mutectJsonLog,
                tumor = pairInfo.tumorId,
                normal = pairInfo.normalId,
                pairName = pairInfo.analysisPairId,
                referenceFa = referenceFa,
                normalFinalBam = pairInfo.normalFinalBam,
                tumorFinalBam = pairInfo.tumorFinalBam,
                highMem = mutect2HighMem,
                mutect2ChrXRawVcfs = mutect2ChrXRawVcfs,
                mutect2ChrXRawStats = mutect2ChrXRawStats
        }
    }
    File mutect2Final = select_first([Mutect2.mutect2, mutect2Preexist])
    if (library == 'WGS') {
        if (!defined(candidateSmallIndelsPreexist)){
            call manta.Manta {
                input:
                    mantaJsonLog = mantaJsonLog,
                    tumor = pairInfo.tumorId,
                    normal = pairInfo.normalId,
                    callRegions = callRegions,
                    referenceFa = referenceFa,
                    pairName = pairInfo.analysisPairId,
                    normalFinalBam = pairInfo.normalFinalBam,
                    tumorFinalBam = pairInfo.tumorFinalBam,
                    highMem = mantaHighMem
            }
        }
        if (!defined(strelka2SnvPreexist)) {
            IndexedVcf candidateSmallIndelsTemp = select_first([Manta.candidateSmallIndels, candidateSmallIndelsPreexist])
            call strelka2.Strelka2 as strelka2Wgs {
                input:
                    library = library,
                    strelkaJsonLog = strelkaJsonLog,
                    configureStrelkaSomaticWorkflow = configureStrelkaSomaticWorkflow,
                    tumor = pairInfo.tumorId,
                    normal = pairInfo.normalId,
                    callRegions = callRegions,
                    intervalListBed = intervalListBed,
                    candidateSmallIndels = candidateSmallIndelsTemp,
                    referenceFa = referenceFa,
                    pairName = pairInfo.analysisPairId,
                    normalFinalBam = pairInfo.normalFinalBam,
                    tumorFinalBam = pairInfo.tumorFinalBam
            }
        }
    }
    
    if (library == 'Exome') {
        call utils.CreateBlankFile as createVcf {
            input:
                fileId = "~{pairInfo.analysisPairId}_nonMantaVcf_"
        }
        
        call utils.CreateBlankFile as createVcfIndex {
            input:
                fileId = "~{pairInfo.analysisPairId}_nonMantaVcfIndex_"
        }
        
        call utils.CreateBlankFile as createfilteredMantaSV {
            input:
                fileId = "~{pairInfo.analysisPairId}_filteredMantaSV_"
        }
        
        IndexedVcf nonMantaCandidateSmallIndels = object {
                vcf : createVcf.blankFile,
                index : createVcfIndex.blankFile
            }
        File blankFilteredMantaSV = createfilteredMantaSV.blankFile
        if (!defined(strelka2SnvPreexist)) {
            call strelka2.Strelka2 as strelka2Exome {
                input:
                    library = library,
                    strelkaJsonLog = strelkaJsonLog,
                    configureStrelkaSomaticWorkflow = configureStrelkaSomaticWorkflow,
                    tumor = pairInfo.tumorId,
                    normal = pairInfo.normalId,
                    callRegions = callRegions,
                    intervalListBed = intervalListBed,
                    candidateSmallIndels = nonMantaCandidateSmallIndels,
                    referenceFa = referenceFa,
                    pairName = pairInfo.analysisPairId,
                    normalFinalBam = pairInfo.normalFinalBam,
                    tumorFinalBam = pairInfo.tumorFinalBam
            }
        }
    }
    
    IndexedVcf candidateSmallIndelsFinal = select_first([Manta.candidateSmallIndels, nonMantaCandidateSmallIndels, candidateSmallIndelsPreexist])
    IndexedVcf strelka2SnvsFinal = select_first([strelka2Wgs.strelka2Snvs, strelka2Exome.strelka2Snvs, strelka2SnvsPreexist])
    IndexedVcf strelka2IndelsFinal = select_first([strelka2Wgs.strelka2Indels, strelka2Exome.strelka2Indels, strelka2IndelsPreexist])
    File strelka2SnvFinal = select_first([strelka2Wgs.strelka2Snv, strelka2Exome.strelka2Snv, strelka2SnvPreexist])
    File strelka2IndelFinal = select_first([strelka2Wgs.strelka2Indel, strelka2Exome.strelka2Indel, strelka2IndelPreexist])
    if (!defined(lancetPreexist)){
        call lancet.Lancet {
            input:
                library = library,
                local = local,
                splitBedsWgsFewNodes = splitBedsWgsFewNodes,
                lancetJsonLog = lancetJsonLog,
                tumorId = pairInfo.tumorId,
                normalId = pairInfo.normalId,
                listOfChroms = listOfChroms,
                splitBedsWgs = splitBedsWgs,
                chromBeds = chromBeds,
                referenceFa = referenceFa,
                pairName = pairInfo.analysisPairId,
                normalFinalBam = pairInfo.normalFinalBam,
                tumorFinalBam = pairInfo.tumorFinalBam,
                serviceAccountKey = serviceAccountKey,
                gcpProject = gcpProject
        }
    }
    File lancetFinal = select_first([Lancet.lancet, lancetPreexist])
    if (!defined(gridssVcfPreexist)){
        call gridss.Gridss {
            input:
                tumor = pairInfo.tumorId,
                normal = pairInfo.normalId,
                pairName = pairInfo.analysisPairId,
                bwaReference = bwaReference,
                gridssAdditionalReference = gridssAdditionalReference,
                listOfChroms = listOfChroms,
                normalFinalBam = pairInfo.normalFinalBam,
                tumorFinalBam = pairInfo.tumorFinalBam,
                bsGenome = bsGenome,
                ponTarGz = ponTarGz,
                highMem = gridssHighMem,
                preMemoryGb = gridssPreMemoryGb,
                filterMemoryGb = gridssFilterMemoryGb,
                gridssUnfilteredVcfChromsPreexist = gridssUnfilteredVcfChromsPreexist
        }
    }
    IndexedVcf gridssVcfFinal = select_first([Gridss.gridssVcf, gridssVcfPreexist])
    if (library == 'WGS') {
        if (!defined(bicseq2Preexist)){
            call bicseq2.BicSeq2 {
                input:
                    tumor = pairInfo.tumorId,
                    normal = pairInfo.normalId,
                    readLength = readLength,
                    coordReadLength = coordReadLength,
                    uniqCoords = uniqCoords,
                    bicseq2ConfigFile = bicseq2ConfigFile,
                    bicseq2SegConfigFile = bicseq2SegConfigFile,
                    chromFastas = chromFastas,
                    listOfChromsFull = listOfChromsFull,
                    pairName = pairInfo.analysisPairId,
                    referenceFa = referenceFa,
                    normalFinalBam = pairInfo.normalFinalBam,
                    tumorFinalBam = pairInfo.tumorFinalBam,
                    tumorMedianInsertSize = tumorMedianInsertSize,
                    normalMedianInsertSize = normalMedianInsertSize,
                    tumorInsertSizeMetrics = tumorInsertSizeMetrics,
                    normalInsertSizeMetrics = normalInsertSizeMetrics,
                    lambda = lambda
            }
        }
        File bicseq2Final = select_first([BicSeq2.bicseq2, bicseq2Preexist])
        File bicseq2PngFinal = select_first([BicSeq2.bicseq2Png, bicseq2PngPreexist])
    }

    output {
        # Mutect2
        File mutect2 = mutect2Final
        # Manta
        IndexedVcf candidateSmallIndels = candidateSmallIndelsFinal
        Array[IndexedVcf] diploidSV = select_all([Manta.diploidSV, diploidSVPreexist])
        Array[IndexedVcf] somaticSV = select_all([Manta.somaticSV, somaticSVPreexist])
        Array[IndexedVcf]  candidateSV = select_all([Manta.candidateSV, candidateSVPreexist])
        Array[File] unfilteredMantaSV = select_all([Manta.unfilteredMantaSV, unfilteredMantaSVPreexist])
        File filteredMantaSV = select_first([Manta.filteredMantaSV, filteredMantaSVPreexist, blankFilteredMantaSV])
        # Strelka2
        IndexedVcf strelka2Snvs = strelka2SnvsFinal
        IndexedVcf strelka2Indels = strelka2IndelsFinal
        File strelka2Snv = strelka2SnvFinal
        File strelka2Indel = strelka2IndelFinal
        # Lancet
        File lancet = lancetFinal
        # Gridss
        IndexedVcf gridssVcf = gridssVcfFinal
        # Bicseq2
        File? bicseq2Png = bicseq2PngFinal
        File? bicseq2 = bicseq2Final
    }
}