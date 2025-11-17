version 1.0

import "alignment_analysis/kourami_wfk.wdl" as kourami
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

workflow AlignmentAnalysis {
    # command
    #   Call variants in BAMs
    #   merge and filter raw VCFs
    #   annotate
    input {
        Array[SampleBamInfo]+ sampleBamInfos
        
        # kourami
        BwaReference kouramiReference
        File kouramiFastaGem1Index
        IndexedReference referenceFa
        
        # needed for work on cloud without localizing
        File? serviceAccountKey
        String? gcpProject

    }
    
    scatter(sampleBamInfo in sampleBamInfos) {
        call kourami.Kourami {
            input:
                sampleId=sampleBamInfo.sampleId,
                kouramiReference=kouramiReference,
                finalBam=sampleBamInfo.finalBam,
                kouramiFastaGem1Index=kouramiFastaGem1Index,
                referenceFa=referenceFa,
                serviceAccountKey = serviceAccountKey,
                gcpProject = gcpProject
        }
    }
    
    output {
        Array[File] kouramiResult = Kourami.result
    }
}
