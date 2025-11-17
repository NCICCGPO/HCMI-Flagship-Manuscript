version 1.0
import "merge_vcf.wdl" as mergeVcf
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
#    Nico Robine (nrobine@nygenome.org)
#    Minita Shah (mshah@nygenome.org)
#    Timothy Chu (tchu@nygenome.org)
#    Will Hooper (whooper@nygenome.org)
#
# ================== /COPYRIGHT ===============================================

workflow MakePairHighConfidenceVcf {
    input {
        File intersectionParticipantVcf
        String tumorBarcodeAliquot
        String normalBarcodeAliquot
        String participantId
        IndexedReference referenceFa
    }
        
    call mergeVcf.CompressIndexVcf {
    input:
        vcf = intersectionParticipantVcf,
        memoryGb = 4
    }
    
    call mergeVcf.SelectVcfSamples {
                input:
                    retainedSamples = [normalBarcodeAliquot, tumorBarcodeAliquot],
                    centerVcf = intersectionParticipantVcf,
                    referenceFa = referenceFa,
                    diskSize = 30
    }
    
    Int diskSize = ceil( size(SelectVcfSamples.nonNormalVcf, "GB") * 2) + ceil(size(intersectionParticipantVcf)) + 5
              
    call mergeVcf.FilterPairHighConfidence {
                input:
                    intersectionPairVcf = SelectVcfSamples.nonNormalVcf,
                    outputVcfPath = "~{tumorBarcodeAliquot}--~{normalBarcodeAliquot}.plusReadSupportAndM2MultiPass.intersection.vcf",
                    intersectionParticipantVcf = CompressIndexVcf.vcfCompressedIndexed,
                    tumorBarcodeAliquot = tumorBarcodeAliquot,
                    normalBarcodeAliquot = normalBarcodeAliquot,
                    participantId = participantId,
                    diskSize = diskSize
            }
    
    output {
        File intersectionPairHighConfidencVcf = FilterPairHighConfidence.intersectionPairHighConfidencVcf 
    }
}
