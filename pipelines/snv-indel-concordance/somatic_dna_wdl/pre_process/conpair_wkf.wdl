version 1.0

import "../wdl_structs.wdl"
import "qc.wdl"

workflow Conpair {
    # command
    input {
        File tumorPileupsConpair
        File normalPileupsConpair
        String tumor
        String normal
        String pairName
        String analysisPairId
        String analysisTumorId

        File markerTxtFile

        Int threads=1
        Int memoryGb=4
        String qcDir = "Sample_~{analysisTumorId}/qc"
    }

    call qc.VerifyConcordanceAll {
        input:
            pileupsTumorConpair = tumorPileupsConpair,
            pileupsNormalConpair = normalPileupsConpair,
            markerTxtFile = markerTxtFile,
            pairName = pairName,
            memoryGb = memoryGb,
            threads = threads,
            concordanceAllPath = "~{qcDir}/~{analysisPairId}.concordance.all.conpair.txt"
    }

    call qc.VerifyConcordanceHomoz {
        input:
            pileupsTumorConpair = tumorPileupsConpair,
            pileupsNormalConpair = normalPileupsConpair,
            markerTxtFile = markerTxtFile,
            pairName = pairName,
            memoryGb = memoryGb,
            threads = threads,
            concordanceHomozPath = "~{qcDir}/~{analysisPairId}.concordance.homoz.conpair.txt"
    }

    call qc.Contamination {
        input:
            pileupsTumorConpair = tumorPileupsConpair,
            pileupsNormalConpair = normalPileupsConpair,
            markerTxtFile = markerTxtFile,
            pairName = pairName,
            memoryGb = memoryGb,
            threads = threads,
            contaminationPath = "~{qcDir}/~{analysisPairId}.contamination.conpair.txt"
    }

    output {
        File concordanceAll = VerifyConcordanceAll.concordanceAll
        File concordanceHomoz = VerifyConcordanceHomoz.concordanceHomoz
        File contamination = Contamination.contamination
    }
}
