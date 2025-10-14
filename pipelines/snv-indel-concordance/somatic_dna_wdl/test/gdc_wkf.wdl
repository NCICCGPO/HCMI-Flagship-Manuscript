version 1.0

import "tests.wdl" as tests
import "../wdl_structs.wdl"
import "../merge_vcf/merge_vcf.wdl" as mergeVcf


workflow GdcComparePair {
    input {
        String name = "gdc_comparison"
        
        String pairId
        String tumor
        String normal
        File nygcVcf
        File gdcSnvVcf
        File gdcIndelVcf
        File nygcFinalBedPe
        File gdcSvVcf
    
        String vepGenomeBuild
        IndexedReference referenceFa
        Array[String] listOfChroms
    }
 
        
    call mergeVcf.Gatk4MergeSortVcf {
        input:
            tempVcfs = [gdcSnvVcf, gdcIndelVcf],
            sortedVcfPath = "~{pairId}.gdc.vcf",
            referenceFa = referenceFa,
            memoryGb = 8,
            diskSize = 10
    }
    
    call tests.FilterHighConfidence {
        input:
            vcf = nygcVcf,
            filteredVcfPath = "~{pairId}.nygc.filtered.vcf",
    }
    
    call tests.SomPy {
        input:
            referenceFa = referenceFa,
            pairId = pairId,
            oldVcf = Gatk4MergeSortVcf.sortedVcf.vcf,
            newVcf = FilterHighConfidence.filteredVcf
    }
    
    call tests.SummarizeVcf as finalSummarizeVcf {
        input:
            pairId = pairId,
            name = "~{name}.vcf",
            oldOnlyVcf = SomPy.oldOnlyVcf,
            newOnlyVcf = SomPy.newOnlyVcf,
            concordantVcf = SomPy.concordantVcf
    }
    
    call tests.VcfToBedPe {
        input:
            tumor = tumor,
            normal = normal,
            vcf = gdcSvVcf,
            outFileBedpePath = "~{pairId}.sv.gdc.bedpe",
            vepGenomeBuild = vepGenomeBuild,
            listOfChroms = listOfChroms
    }
    
    Boolean runSv = false
    if (runSv) {
    
        call tests.CompareBedPe as finalCompareBedPe {
            input:
                pairId = pairId,
                oldBedpe = VcfToBedPe.outFileBedpe,
                newBedpe = nygcFinalBedPe
        }
        
        call tests.CompareSvGenes as finalCompareSvGenes {
            input:
                pairId = pairId,
                name = name,
                concordanceBedPe = finalCompareBedPe.outFileBedpe
        }
    }

    output {
        #File finalBedTable = finalCompareBedPe.outFileSummary
        #File finalBedpeJoined = finalCompareBedPe.outFileBedpe
        #File svGenesTable = finalCompareSvGenes.svGenesTable
        File summaryOutputTable = finalSummarizeVcf.summaryOutputTable
        File detailedOutputTable = finalSummarizeVcf.detailedOutputTable

    }
}

