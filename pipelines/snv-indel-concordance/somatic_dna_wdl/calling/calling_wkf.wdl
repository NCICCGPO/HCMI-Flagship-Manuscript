version 1.0

import "mutect2_wkf.wdl" as mutect2
import "strelka2_wkf.wdl" as strelka2
import "manta_wkf.wdl" as manta
import "lancet_wkf.wdl" as lancet
import "gridss_wkf.wdl" as gridss
import "bicseq2_wkf.wdl" as bicseq2
import "../tasks/utils.wdl" as utils

import "../wdl_structs.wdl"

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



workflow Calling {
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
        Int tumorMedianInsertSize = 400
        Int normalMedianInsertSize = 400
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

    }
    call mutect2.Mutect2 {
        input:
            local = local,
            library = library,
            invertedIntervalListBed = invertedIntervalListBed,
            callerIntervals = callerIntervals,
            mutectJsonLogFilter = mutectJsonLogFilter,
            mutectJsonLog = mutectJsonLog,
            tumor = pairInfo.tumorId,
            normal = pairInfo.normalId,
            pairName = pairInfo.pairId,
            referenceFa = referenceFa,
            normalFinalBam = pairInfo.normalFinalBam,
            tumorFinalBam = pairInfo.tumorFinalBam,
            highMem = mutect2HighMem,
            mutect2ChrXRawVcfs = mutect2ChrXRawVcfs,
            mutect2ChrXRawStats = mutect2ChrXRawStats
    }

    if (library == 'WGS') {
        call manta.Manta {
            input:
                mantaJsonLog = mantaJsonLog,
                tumor = pairInfo.tumorId,
                normal = pairInfo.normalId,
                callRegions = callRegions,
                referenceFa = referenceFa,
                pairName = pairInfo.pairId,
                normalFinalBam = pairInfo.normalFinalBam,
                tumorFinalBam = pairInfo.tumorFinalBam,
                highMem = mantaHighMem
        }
        
        call strelka2.Strelka2 as strelka2Wgs {
            input:
                library = library,
                strelkaJsonLog = strelkaJsonLog,
                configureStrelkaSomaticWorkflow = configureStrelkaSomaticWorkflow,
                tumor = pairInfo.tumorId,
                normal = pairInfo.normalId,
                callRegions = callRegions,
                intervalListBed = intervalListBed,
                candidateSmallIndels = Manta.candidateSmallIndels,
                referenceFa = referenceFa,
                pairName = pairInfo.pairId,
                normalFinalBam = pairInfo.normalFinalBam,
                tumorFinalBam = pairInfo.tumorFinalBam
        }
    }
    
    if (library == 'Exome') {
        call utils.CreateBlankFile as createVcf {
            input:
                fileId = "~{pairInfo.pairId}_nonMantaVcf_"
        }
        
        call utils.CreateBlankFile as createVcfIndex {
            input:
                fileId = "~{pairInfo.pairId}_nonMantaVcfIndex_"
        }
        
        IndexedVcf nonMantaCandidateSmallIndels = object {
                vcf : createVcf.blankFile,
                index : createVcfIndex.blankFile
            }
            
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
                pairName = pairInfo.pairId,
                normalFinalBam = pairInfo.normalFinalBam,
                tumorFinalBam = pairInfo.tumorFinalBam
        }
    }
    
    IndexedVcf candidateSmallIndelsFinal = select_first([Manta.candidateSmallIndels, nonMantaCandidateSmallIndels])
    IndexedVcf strelka2SnvsFinal = select_first([strelka2Wgs.strelka2Snvs, strelka2Exome.strelka2Snvs])
    IndexedVcf strelka2IndelsFinal = select_first([strelka2Wgs.strelka2Indels, strelka2Exome.strelka2Indels])
    File strelka2SnvFinal = select_first([strelka2Wgs.strelka2Snv, strelka2Exome.strelka2Snv])
    File strelka2IndelFinal = select_first([strelka2Wgs.strelka2Indel, strelka2Exome.strelka2Indel])
    

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
            pairName = pairInfo.pairId,
            normalFinalBam = pairInfo.normalFinalBam,
            tumorFinalBam = pairInfo.tumorFinalBam,
            serviceAccountKey = serviceAccountKey,
            gcpProject = gcpProject
    }

    call gridss.Gridss {
        input:
            tumor = pairInfo.tumorId,
            normal = pairInfo.normalId,
            pairName = pairInfo.pairId,
            bwaReference = bwaReference,
            gridssAdditionalReference = gridssAdditionalReference,
            listOfChroms = listOfChroms,
            normalFinalBam = pairInfo.normalFinalBam,
            tumorFinalBam = pairInfo.tumorFinalBam,
            bsGenome = bsGenome,
            ponTarGz = ponTarGz,
            highMem = gridssHighMem,
            preMemoryGb = gridssPreMemoryGb,
            filterMemoryGb = gridssFilterMemoryGb
    }
    
    if (library == 'WGS') {

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
                pairName = pairInfo.pairId,
                referenceFa = referenceFa,
                normalFinalBam = pairInfo.normalFinalBam,
                tumorFinalBam = pairInfo.tumorFinalBam,
                tumorMedianInsertSize = tumorMedianInsertSize,
                normalMedianInsertSize = normalMedianInsertSize,
                lambda = lambda
        }
    }

    output {
        # Mutect2
        File mutect2 = Mutect2.mutect2
        # Manta
        IndexedVcf candidateSmallIndels = candidateSmallIndelsFinal
        IndexedVcf? diploidSV = Manta.diploidSV
        IndexedVcf?  somaticSV = Manta.somaticSV
        IndexedVcf?  candidateSV = Manta.candidateSV
        File? unfilteredMantaSV = Manta.unfilteredMantaSV
        File? filteredMantaSV = Manta.filteredMantaSV
        # Strelka2
        IndexedVcf strelka2Snvs = strelka2SnvsFinal
        IndexedVcf strelka2Indels = strelka2IndelsFinal
        File strelka2Snv = strelka2SnvFinal
        File strelka2Indel = strelka2IndelFinal
        # Lancet
        File lancet = Lancet.lancet
        # Gridss
        IndexedVcf gridssVcf = Gridss.gridssVcf
        # Bicseq2
        File? bicseq2Png = BicSeq2.bicseq2Png
        File? bicseq2 = BicSeq2.bicseq2
    }
}
