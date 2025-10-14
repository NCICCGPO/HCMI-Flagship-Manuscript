version 1.0

import "merge_vcf/merge_vcf_exome_wgs_wkf.wdl" as mergeVcf
import "wdl_structs.wdl"

# ================== COPYRIGHT ================================================
# New York Genome Center
# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright (2023) by the New York
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

workflow MergeVcf {
    # command
    #   Call variants in BAMs
    #   merge and filter raw VCFs
    #   annotate
    input {
        Boolean external = false
        File wgsVcf
        File exomeVcf
        String tumorId
        String normalId
        String pairId
        String exomeTumorId
        String exomeNormalId
        String exomePairId
        
        IndexedReference referenceFa
        Array[String]+ listOfChroms
        
    }
    
    call mergeVcf.MergeVcfExomeWgs {
        input:
            external = external,
            wgsVcf = wgsVcf,
            exomeVcf = exomeVcf,
            tumorId=tumorId,
            normalId=normalId,
            pairId=pairId,
            exomeTumorId=exomeTumorId,
            exomeNormalId=exomeNormalId,
            exomePairId=exomePairId,
            referenceFa = referenceFa,
            listOfChroms = listOfChroms
    }
    
    output {
        File mergedVcf = MergeVcfExomeWgs.mergedVcf
    }
}
