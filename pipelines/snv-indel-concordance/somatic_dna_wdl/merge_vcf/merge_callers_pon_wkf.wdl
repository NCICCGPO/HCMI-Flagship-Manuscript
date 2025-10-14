version 1.0

import "merge_vcf.wdl" as mergeVcf
import "../wdl_structs.wdl"

# note that we will still need to match dictionaries or provide matching reference files

workflow MergeCallersPon {
    input {
        String tumor
        Array[String]+ listOfChroms
        Array[IndexedVcf]+ allVcfCompressed
    }

    scatter(vcfCompressed in allVcfCompressed) {
            File allVcfCompressedFile = vcfCompressed.vcf
        }
    Array[File]+ allVcfCompressedList = allVcfCompressedFile

    scatter(chrom in listOfChroms) {
        call mergeVcf.MergeCallers as allMergeCallers {
            input:
                chrom = chrom,
                pairName = tumor,
                allVcfCompressed = allVcfCompressed,
                allVcfCompressedList = allVcfCompressedList,
        }

        #  =================================================================
        #                     Merge columns
        #  =================================================================

        call mergeVcf.MergeColumnsPon {
            input:
                chrom = chrom,
                tumor = tumor,
                supportedChromVcf = allMergeCallers.mergedChromVcf
        }
    }

    output {
            Array[File] finalChromVcf = MergeColumnsPon.columnChromVcf
        }
}
