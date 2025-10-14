version 1.0

import "alignment_analysis.wdl" as alignment_analysis
import "../wdl_structs.wdl"

workflow Msi {
    input {
        String normal
        String pairName
        String mantisBedByIntervalListPath = "mantisBedByIntervalListPath.bed"
        File mantisBed
        File intervalListBed
        IndexedReference referenceFa
        Bam normalFinalBam
        Bam tumorFinalBam
    }
    
    call alignment_analysis.BedtoolsIntersect {
        input:
            mantisBedByIntervalListPath = mantisBedByIntervalListPath,
            mantisBed = mantisBed,
            intervalListBed = intervalListBed
    }
    
    call alignment_analysis.MantisExome {
        input:
            tumorFinalBam = tumorFinalBam,
            normalFinalBam = normalFinalBam,
            pairName = pairName,
            mantisBedByIntervalList = BedtoolsIntersect.mantisBedByIntervalList,
            referenceFa = referenceFa
    }
    
    call alignment_analysis.MantisRethreshold {
        input:
            normal = normal,
            pairName = pairName,
            mantisWxsStatus = MantisExome.mantisWxsStatus
    }
    
    output {
        File mantisWxsKmerCountsFinal = MantisExome.mantisWxsKmerCountsFinal
        File mantisWxsKmerCountsFiltered = MantisExome.mantisWxsKmerCountsFiltered
        File mantisExomeTxt = MantisExome.mantisExomeTxt
        File mantisStatusFinal = MantisRethreshold.mantisStatusFinal
    }
}
