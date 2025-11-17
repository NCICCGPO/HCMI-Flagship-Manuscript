version 1.0

import "align_fastq.wdl" as alignFastq
import "../wdl_structs.wdl"

workflow AlignFastq {
    # command
    #   align FASTQ files
    input {
        Array[Fastqs] listOfFastqPairs
        BwaMem2Reference bwamem2Reference
        File adaptersFa
        Boolean trim = true
        # resources
        Int bwamem2Mem
        Int threads
        Int bwamem2Threads
        Int additionalDiskSize = 20
    }

    scatter(fastqs in listOfFastqPairs) {
        Int fastqsSize = ceil(size(fastqs.fastqR1, "GB") + size (fastqs.fastqR2, "GB"))

        # Assumption: The bam will be about the same size as the input fastqs. So double
        # the size to account for input and output.
        Int skewerDiskSize = (2 * fastqsSize) + additionalDiskSize
        Int alignDiskSize = (3 * fastqsSize) + additionalDiskSize

        if (trim) {
            call alignFastq.Skewer {
                input:
                    fastqs = fastqs,
                    adaptersFa = adaptersFa,
                    diskSize = skewerDiskSize
            }
        }

        Fastqs fastqsAlign = select_first([Skewer.skewerFastqs, fastqs])

        call alignFastq.AlignBwaMem2 {
            input:
                fastqsAlign = fastqsAlign,
                bwamem2Reference = bwamem2Reference,
                memoryGb = bwamem2Mem,
                threads = threads,
                bwamem2Threads = bwamem2Threads,
                diskSize = alignDiskSize
        }

        Int shortMarkDiskSize = ceil( 4 *  size(AlignBwaMem2.laneBam, "GB")) + additionalDiskSize
        call alignFastq.ShortAlignMark {
            input:
                laneBam = AlignBwaMem2.laneBam,
                bamBase = fastqs.readGroupId,
                memoryGb = 16,
                diskSize = shortMarkDiskSize
        }

        Int fixmateDiskSize = ceil( 4 *  size(ShortAlignMark.laneBamMark, "GB")) + additionalDiskSize
        call alignFastq.Fixmate {
            input:
                laneBamMark = ShortAlignMark.laneBamMark,
                bamBase = fastqs.readGroupId,
                memoryGb = 8,
                diskSize = fixmateDiskSize
        }
    }

    output {
        Array[File] laneFixmateBam = Fixmate.laneFixmateBam
    }
}
