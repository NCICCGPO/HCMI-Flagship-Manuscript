version 1.0

import "mutect2_wkf.wdl" as mutect2
import "strelka2_wkf.wdl" as strelka2
import "manta_wkf.wdl" as manta
import "lancet_wkf.wdl" as lancet
import "../tasks/utils.wdl" as utils

import "../wdl_structs.wdl"

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
        File intervalListBed
        File invertedIntervalListBed
        IndexedReference referenceFa
        #   Manta
        IndexedTable callRegions
        BwaReference bwaReference
        #   Lancet
        Array[File] splitBedsWgs
        Map[String, File] chromBeds

        File lancetJsonLog
        File mantaJsonLog
        File strelkaJsonLog
        File mutectJsonLog
        File mutectJsonLogFilter
        File configureStrelkaSomaticWorkflow

        # Resources need a lot of fine grained control
        Boolean mantaHighMem = false
        Boolean mutect2HighMem = false
        # need to create docker with gcloud
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
            highMem = mutect2HighMem
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
    }
}
