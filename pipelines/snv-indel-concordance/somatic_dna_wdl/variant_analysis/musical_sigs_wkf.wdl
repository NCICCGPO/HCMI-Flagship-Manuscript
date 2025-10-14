version 1.0

import "variant_analysis.wdl" as variant_analysis
import "../test/tests.wdl" as tests
import "../wdl_structs.wdl"

workflow MusicalSig {
    # command
    #   run MusicalSig on filtered VCF
    input {
        String pairId
        File mainVcf
        String vepGenomeBuild
        File cosmicSbsReference
        File cosmicDbsReference
        File cosmicIdReference
    }
    call variant_analysis.SigProfilerMatrixGenerator {
        input:
            filteredVcf = mainVcf,
            pairId = pairId,
            vepGenomeBuild = vepGenomeBuild
    }
    
    call variant_analysis.MusicalRefit {
        input:
            pairId = pairId,
            outputSbs = SigProfilerMatrixGenerator.outputSbs,
            outputDbs = SigProfilerMatrixGenerator.outputDbs,
            outputId = SigProfilerMatrixGenerator.outputId,
            cosmicSbsReference = cosmicSbsReference,
            cosmicDbsReference = cosmicDbsReference,
            cosmicIdReference = cosmicIdReference
    }

    output {
       File outputSbsMatrix = SigProfilerMatrixGenerator.outputSbs
       File outputDbsMatrix = SigProfilerMatrixGenerator.outputDbs
       File outputIdMatrix = SigProfilerMatrixGenerator.outputId
       File sigsSbs = MusicalRefit.sigsSbs
       File? sigsDbs = MusicalRefit.sigsDbs
       File? sigsId = MusicalRefit.sigsId
    }
}