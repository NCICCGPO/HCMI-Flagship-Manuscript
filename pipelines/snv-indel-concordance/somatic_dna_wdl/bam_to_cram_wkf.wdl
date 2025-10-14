version 1.0

import "./wdl_structs.wdl"
import "tasks/bam_cram_conversion.wdl" as cramConversion


workflow BamToCramWorkflow {
    input {
        IndexedReference referenceFa
        Array[PairInfo]+ pairInfos
        Array[SampleBamInfo]+ normalSampleBamInfos
    }

    call cramConversion.UniqueBams as uniqueBams {
        input:
            pairInfosJson = write_json(pairInfos)
    }

    scatter(bamInfo in uniqueBams.uniqueBams) {
        call cramConversion.SamtoolsBamToCram as bamToCram {
            input:
                inputBam = bamInfo.finalBam,
                referenceFa = referenceFa,
                sampleId = bamInfo.sampleId,
                diskSize = (ceil(size(bamInfo.finalBam.bam, "GB") * 1.5)) + 50
        }
    }

    output {
        Array[SampleCramInfo] crams = bamToCram.cramInfo

    }
}
