version 1.0
import "../wdl_structs.wdl"
import "../tasks/utils.wdl" as utils
import "../pre_process/qc.wdl" as qc
import "../pre_process/conpair_wkf.wdl" as conpair
import "../baf/baf_wkf.wdl" as baf
import "../alignment_analysis/msi_wkf.wdl" as msi
import "../alignment_analysis/telomeasures_wkf.wdl" as telomeasures

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


workflow NygcPairedTasks {
    input {
        # mantis
        File mantisBed
        File intervalListBed
        IndexedReference referenceFa
        File markerTxtFile
        File markerBedFile
        File teloTargetIndexGem
        Float maxMismatches = 0.2
        File adaptersFa
        Boolean trim = false
        # TelFusDetector parameters
        Int purity = 1
        # patient-specific inputs
        PairInfo pairInfo
        # need to create docker with gcloud
        File? serviceAccountKey
        String? gcpProject
        File? conpairPileupNormalPreexist
        File? conpairPileupTumorPreexist
    }
    if (!defined(conpairPileupNormalPreexist)) {
        String normalPileupsConpairPath = "~{pairInfo.analysisNormalId}_pileups_table.txt"
        call qc.ConpairPileup as conpairPileupNormalRun {
                input:
                    markerBedFile = markerBedFile,
                    referenceFa = referenceFa,
                    finalBam = pairInfo.normalFinalBam,
                    sampleId = pairInfo.normalId,
                    pileupsConpairPath = normalPileupsConpairPath,
                    memoryGb = 4,
                    threads = 1,
                    diskSize = ceil(size(pairInfo.normalFinalBam.bam, "GB"))  + 20
        }
    }
    File conpairPileupNormalFinal = select_first([conpairPileupNormalRun.pileupsConpair, conpairPileupNormalPreexist])
    if (!defined(conpairPileupTumorPreexist)) {
        String tumorPileupsConpairPath = "~{pairInfo.analysisTumorId}_pileups_table.txt"
        call qc.ConpairPileup as conpairPileupTumorRun {
            input:
                markerBedFile = markerBedFile,
                referenceFa = referenceFa,
                finalBam = pairInfo.tumorFinalBam,
                sampleId = pairInfo.tumorId,
                pileupsConpairPath = tumorPileupsConpairPath,
                memoryGb = 4,
                threads = 1,
                diskSize = ceil(size(pairInfo.tumorFinalBam.bam, "GB"))  + 20
        }
    }
    File conpairPileupTumorFinal = select_first([conpairPileupTumorRun.pileupsConpair, conpairPileupTumorPreexist])
    call conpair.Conpair {
                input:
                    tumorPileupsConpair = conpairPileupTumorFinal,
                    normalPileupsConpair = conpairPileupNormalFinal,
                    analysisPairId = pairInfo.analysisPairId,
                    analysisTumorId = pairInfo.analysisTumorId,
                    tumor = pairInfo.tumorId,
                    normal = pairInfo.normalId,
                    pairName = pairInfo.pairId,
                    markerTxtFile = markerTxtFile
    }
    call telomeasures.Telomeasures as tumorTelomeasures {
            input:
                finalBam = pairInfo.tumorFinalBam,
                referenceFa = referenceFa,
                sampleId = pairInfo.tumorId,
                analysisId = pairInfo.analysisTumorId,
                teloTargetIndexGem = teloTargetIndexGem,
                maxMismatches = maxMismatches,
                trim = trim,
                adaptersFa = adaptersFa,
                purity = purity,
                serviceAccountKey = serviceAccountKey,
                gcpProject = gcpProject
    }
    call telomeasures.Telomeasures as normalTelomeasures {
            input:
                finalBam = pairInfo.normalFinalBam,
                referenceFa = referenceFa,
                sampleId = pairInfo.normalId,
                analysisId = pairInfo.analysisNormalId,
                teloTargetIndexGem = teloTargetIndexGem,
                maxMismatches = maxMismatches,
                trim = trim,
                adaptersFa = adaptersFa,
                purity = 1,
                serviceAccountKey = serviceAccountKey,
                gcpProject = gcpProject
    }
    call msi.Msi {
                input:
                    normal = pairInfo.normalId,
                    pairName = pairInfo.analysisPairId,
                    mantisBed = mantisBed,
                    intervalListBed = intervalListBed,
                    referenceFa = referenceFa,
                    tumorFinalBam = pairInfo.tumorFinalBam,
                    normalFinalBam = pairInfo.normalFinalBam
            }
    output {
        # Conpair
        File concordanceAll = Conpair.concordanceAll
        File concordanceHomoz = Conpair.concordanceHomoz
        File contamination = Conpair.contamination
        File conpairPileupTumor = conpairPileupTumorFinal
        File conpairPileupNormal = conpairPileupNormalFinal
        # TeloRead Descriptions
        File normalGemLength = normalTelomeasures.gemLength
        # Aligned/Unaligned Telomeric reads
        File normalTeloR1Fastq = normalTelomeasures.teloR1Fastq
        File normalTeloR2Fastq = normalTelomeasures.teloR2Fastq
        File normalSingMappedFastq = normalTelomeasures.singMappedFastq
        File normalAlignmentCsv = normalTelomeasures.alignmentCsv
        # normal teloreads
        # TelFusDetector
        File normalFusionsFiltered = normalTelomeasures.fusionsFiltered
        File normalFusionsPass = normalTelomeasures.fusionsPass
        File normalAllChromosomesCov = normalTelomeasures.allChromosomesCov
        File normalFusionsRates = normalTelomeasures.fusionsRates
        # telomeasures
        File telomeasuresNormal = normalTelomeasures.telomeasures
        # tumor teloreads
        # TeloRead Descriptions
        File tumorGemLength = tumorTelomeasures.gemLength
        # Aligned/Unaligned Telomeric reads
        File tumorTeloR1Fastq = tumorTelomeasures.teloR1Fastq
        File tumorTeloR2Fastq = tumorTelomeasures.teloR2Fastq
        File tumorSingMappedFastq = tumorTelomeasures.singMappedFastq
        File tumorAlignmentCsv = tumorTelomeasures.alignmentCsv
        # TelFusDetector
        File tumorFusionsFiltered = tumorTelomeasures.fusionsFiltered
        File tumorFusionsPass = tumorTelomeasures.fusionsPass
        File tumorAllChromosomesCov = normalTelomeasures.allChromosomesCov
        File tumorFusionsRates = tumorTelomeasures.fusionsRates
        # telomeasures
        File telomeasuresTumor = tumorTelomeasures.telomeasures
        # MSI
        File mantisWxsKmerCountsFinal = Msi.mantisWxsKmerCountsFinal
        File mantisWxsKmerCountsFiltered = Msi.mantisWxsKmerCountsFiltered
        File mantisExomeTxt = Msi.mantisExomeTxt
        File mantisStatusFinal = Msi.mantisStatusFinal
    }
}
