version 1.0

import "test/tests.wdl" as tests
import "wdl_structs.wdl"
import "test/gdc_wkf.wdl" as gdcComparePair


struct GdcRelationship {
    String pairId
    String tumor
    String normal
    File nygcVcf
    File gdcSnvVcf
    File gdcIndelVcf
    File gdcSvVcf
    File nygcFinalBedPe
}


workflow GdcCompare {
    input {
        String name = "gdc_comparison"
        File chromLengths
        File cosmicCensus
        Array[String] listOfChroms
        
        String vepGenomeBuild
        Array[GdcRelationship] gdcRelationships
        IndexedReference referenceFa
    }
 
    
    scatter (gdcRelationship in gdcRelationships) {
        
        call gdcComparePair.GdcComparePair {
            input:
                pairId = gdcRelationship.pairId,
                tumor = gdcRelationship.tumor,
                normal = gdcRelationship.normal,
                nygcVcf = gdcRelationship.nygcVcf,
                gdcIndelVcf = gdcRelationship.gdcIndelVcf,
                gdcSnvVcf = gdcRelationship.gdcSnvVcf,
                gdcSvVcf = gdcRelationship.gdcSvVcf,
                nygcFinalBedPe = gdcRelationship.nygcFinalBedPe,
                listOfChroms = listOfChroms,
                referenceFa = referenceFa,
                vepGenomeBuild = vepGenomeBuild
        }
    }
    
    
    call tests.ConcateTables as finalVcfConcateTables {
        input:
            tables = GdcComparePair.summaryOutputTable,
            outputTablePath = "~{name}.finalVcf.csv"
    }
    
    #Boolean runSv = false
    #if (runSv) {
    #    call tests.ConcateTables as finalBedPeConcateTables {
    #        input:
    #            tables = GdcComparePair.finalBedTable,
    #            outputTablePath = "~{name}.finalBedPe.tsv"
    #    }
    #    
    #    call tests.ConcateTables as finalSvGenesConcateTables {
    #        input:
    #            tables = GdcComparePair.svGenesTable,
    #            outputTablePath = "~{name}.finalSvGenes.csv"
    #    }
    #}
        
    output {
        #File finalBedPeTable = finalBedPeConcateTables.outputTable
        #Array[File] finalBedpeJoined = GdcComparePair.finalBedpeJoined
        #File finalSvGenesTable = finalSvGenesConcateTables.outputTable
        File finalVcfTable = finalVcfConcateTables.outputTable
        Array[File] detailedOutputVcfTable = GdcComparePair.detailedOutputTable

    }
}

