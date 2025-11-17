version 1.0

import "merge_vcf.wdl" as mergeVcf
import "../calling/calling.wdl" as calling
import "../wdl_structs.wdl"

# note that we will still need to match dictionaries or provide matching reference files

workflow MergeCallers {
    input {
        Boolean external = false
        String tumorId
        String normalId
        String pairName
        Array[String]+ listOfChroms
        Array[IndexedVcf]+ allVcfCompressed

        IndexedReference referenceFa
    }

    scatter(vcfCompressed in allVcfCompressed) {
            File allVcfCompressedFile = vcfCompressed.vcf
        }
    Array[File]+ allVcfCompressedList = allVcfCompressedFile

    scatter(chrom in listOfChroms) {
        call mergeVcf.MergeExomeWgsCallersFull {
            input:
                chrom = chrom,
                pairName = pairName,
                tumorId=tumorId,
                normalId=normalId,
                allVcfCompressed = allVcfCompressed,
                allVcfCompressedList = allVcfCompressedList
        }

    }

    output {
            Array[File] mergedChromVcf = MergeExomeWgsCallersFull.columnChromVcf
            Array[File] finalChromVcf = MergeExomeWgsCallersFull.finalChromVcf
        }
}
