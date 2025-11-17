version 1.0

import "alignment_analysis.wdl" as alignmentAnalysis
import "../wdl_structs.wdl"

workflow FastNgsAdmix {
    input {
        Bam normalFinalBam
        File fastNgsAdmixSites
        File fastNgsAdmixSitesBin
        File fastNgsAdmixSitesIdx
        File fastNgsAdmixChroms
        File fastNgsAdmixRef
        File fastNgsAdmixNind
        Int angsdThreads = 4 
        String outprefix
    }

    call alignmentAnalysis.Angsd {
        input:
            normalFinalBam = normalFinalBam,
            fastNgsAdmixSites = fastNgsAdmixSites,
            fastNgsAdmixSitesBin = fastNgsAdmixSitesBin,
            fastNgsAdmixSitesIdx = fastNgsAdmixSitesIdx,
            fastNgsAdmixChroms = fastNgsAdmixChroms,
            threads = angsdThreads,
            outprefix = outprefix
    }

    call alignmentAnalysis.FastNgsAdmix {
        input:
            beagleFile = Angsd.beagleFile,
            fastNgsAdmixRef = fastNgsAdmixRef,
            fastNgsAdmixNind = fastNgsAdmixNind,
            outprefix = outprefix
    }

    output {
        File beagleFile = Angsd.beagleFile
        File fastNgsAdmixQopt = FastNgsAdmix.fastNgsAdmixQopt
    }
}