version 1.0

import "variant_analysis.wdl" as variant_analysis
import "../wdl_structs.wdl"

workflow DeconstructSig {
    # command 
    #   run DeconstructSig on merged VCF
    input {
        String pairId
        File mainVcf
        String vepGenomeBuild
        File cosmicSigs
    }
    
    call variant_analysis.Deconstructsig {
        input:
            mainVcf = mainVcf,
            pairId = pairId,
            cosmicSigs = cosmicSigs,
            vepGenomeBuild = vepGenomeBuild
    }
    
    output {
         File sigs = Deconstructsig.sigs
         File counts = Deconstructsig.counts
         File sigInput = Deconstructsig.sigInput
         File reconstructed = Deconstructsig.reconstructed
         File diff = Deconstructsig.diff
    }
}
    
