version 1.0

import "alignment_analysis.wdl" as telomeasureTasks
import "../wdl_structs.wdl"

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
#
# ================== /COPYRIGHT ===============================================

workflow Telomeasures {
    input {
        Bam finalBam
        String sampleId
        String analysisId
        File teloTargetIndexGem
        IndexedReference referenceFa
        Float maxMismatches = 0.2
        File adaptersFa
        Boolean trim = false
        # TelFusDetector parameters
        Int purity = 1
        Array[String] groupings = ["Chromosomal_region", "Orientation", "Filter"]
        # needed for work on cloud without localizing
        File? serviceAccountKey
        String? gcpProject
        
        Int diskSize = ceil(size(finalBam.bam, "GB")) + 1
    }
    
    # Identify possible TeloReads
    if (trim) {
        call telomeasureTasks.UnclippedPercentTrim {
            input:
                finalBam = finalBam,
                referenceFa = referenceFa,
                sampleId = analysisId,
                teloTargetIndexGem = teloTargetIndexGem,
                maxMismatches = maxMismatches,
                diskSize = diskSize,
                adaptersFa = adaptersFa,
                serviceAccountKey = serviceAccountKey,
                gcpProject = gcpProject
        }
    }
    if (!trim) {
        call telomeasureTasks.UnclippedPercent {
            input:
                finalBam = finalBam,
                referenceFa = referenceFa,
                sampleId = analysisId,
                teloTargetIndexGem = teloTargetIndexGem,
                maxMismatches = maxMismatches,
                diskSize = diskSize,
                serviceAccountKey = serviceAccountKey,
                gcpProject = gcpProject
        }
    }
    File r1MappedFastqRun = select_first([UnclippedPercentTrim.r1MappedFastq, UnclippedPercent.r1MappedFastq])
    File r2MappedFastqRun = select_first([UnclippedPercentTrim.r2MappedFastq, UnclippedPercent.r2MappedFastq])
    File singMappedFastqRun = select_first([UnclippedPercentTrim.singMappedFastq, UnclippedPercent.singMappedFastq])
    File gemLengthRun = select_first([UnclippedPercentTrim.gemLength, UnclippedPercent.gemLength])
    File alignmentCsvRun = select_first([UnclippedPercentTrim.alignmentCsv, UnclippedPercent.alignmentCsv])

    call telomeasureTasks.TelFusDetectorCaller {
        input:
            finalBam = finalBam,
            sampleId = sampleId,
            analysisId = analysisId,
            serviceAccountKey = serviceAccountKey,
            gcpProject = gcpProject
            
    }
    
    call telomeasureTasks.TelFusDetectorRates {
        input:
            fusionsPass = TelFusDetectorCaller.fusionsPass,
            sampleId = sampleId,
            analysisId = analysisId,
            purity = purity,
            groupings = groupings
            
    }
    String telomeasuresPath = "~{analysisId}.telomeasures.summary.csv"
    call telomeasureTasks.TeloRatio {
        input:
            sampleId = sampleId,
            telomeasuresPath = telomeasuresPath,
            alignmentCsv = alignmentCsvRun
    }
    
    output {
        # TeloRead Descriptions
        File gemLength = gemLengthRun
        # Aligned/Unaligned Telomeric reads
        File teloR1Fastq = r1MappedFastqRun
        File teloR2Fastq = r2MappedFastqRun
        File singMappedFastq = singMappedFastqRun
        File alignmentCsv = alignmentCsvRun
        # TelFusDetector
        File fusionsFiltered = TelFusDetectorCaller.fusionsFiltered
        File fusionsPass = TelFusDetectorCaller.fusionsPass
        File allChromosomesCov = TelFusDetectorCaller.allChromosomesCov
        File fusionsRates = TelFusDetectorRates.fusionsRates
        # telomeasures
        File telomeasures = TeloRatio.telomeasures
    }
}