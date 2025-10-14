version 1.0

import "merge_vcf/make_paired_high_confidence_vcf_wkf.wdl" as MakePairHighConfidenceVcf
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

workflow MakePairHighConfidenceVcfs {
    # command
    #   Call variants in BAMs
    #   merge and filter raw VCFs
    #   annotate
    input {
        File intersectionParticipantVcf
        Array[String] tumorBarcodeAliquots
        String normalBarcodeAliquot
        String participantId
        IndexedReference referenceFa
    }
    
    scatter (i in range(length(tumorBarcodeAliquots))) {
        call MakePairHighConfidenceVcf.MakePairHighConfidenceVcf {
            input:
                intersectionParticipantVcf = intersectionParticipantVcf,
                tumorBarcodeAliquot = tumorBarcodeAliquots[i],
                normalBarcodeAliquot = normalBarcodeAliquot,
                participantId = participantId,
                referenceFa = referenceFa
        }
    }
    
    output {
        Array[File] intersectionPairHighConfidencVcf = MakePairHighConfidenceVcf.intersectionPairHighConfidencVcf
    }
}
