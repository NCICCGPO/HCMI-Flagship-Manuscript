version 1.0

import "tests.wdl" as tests
import "../wdl_structs.wdl"


workflow SampleReport {
    input {
        String pipeline = "v7"
        String pairId
        String normal
        String docLink
        File chromLengths
        File cosmicCensus
        IndexedReference referenceFa
        File karyotype
        File navTemplate
        File pandocTemplate
        Array[String] listOfChroms
        
        # All calls
        FinalVcfPairInfo finalVcfPairInfo
        # BAF
        File alleleCountsTxt
        # HC
        File filteredHaplotypecallerAnnotatedVcf
        # MSI
        File mantisStatusFinal
        # Hla
        File kouramiResult
        # Sigs
        File diff
        File sig_input
        File reconstructed
        File sigs
        
        Int bedDiskSize = ceil( size(finalVcfPairInfo.cnvAnnotatedSupplementalBed, "GB") * 2) + 5
        Int bedPeDiskSize = ceil( size(finalVcfPairInfo.svFinalBedPe, "GB") * 2) + 5
        Int vcfDiskSize = ceil( size(finalVcfPairInfo.mainVcf, "GB") * 2) + 5
        Int germVcfDiskSize = ceil( size(filteredHaplotypecallerAnnotatedVcf, "GB") * 2) + 5
        Int bafDiskSize = ceil( size(alleleCountsTxt, "GB") * 2) + 5
    }
    
    call tests.DescribeBed {
        input:
            pairId = finalVcfPairInfo.pairId,
            bed = finalVcfPairInfo.cnvAnnotatedSupplementalBed,
            chromLengths = chromLengths,
            listOfChroms = listOfChroms,
            diskSize = bedDiskSize
    }
    
    call tests.DescribeBedGenes {
        input:
            pairId = finalVcfPairInfo.pairId,
            bed = finalVcfPairInfo.cnvAnnotatedSupplementalBed,
            diskSize = bedDiskSize
    }
    
    call tests.SummarizeFinalVcf {
            input:
                pairId = pairId,
                vcf = finalVcfPairInfo.mainVcf,
                diskSize = vcfDiskSize
        }
        
    call tests.SummarizeFinalGermVcf {
            input:
                sampleId = normal,
                vcf = filteredHaplotypecallerAnnotatedVcf,
                diskSize = germVcfDiskSize
        }
        
    call tests.DescribeBedPeGenes {
        input:
            pairId = finalVcfPairInfo.pairId,
            bedpe = finalVcfPairInfo.svHighConfidenceFinalBedPe,
            diskSize = bedPeDiskSize
    }
    
    call tests.DescribeBedPe as highConfidenceDescribeBedPe {
        input:
            pairId = finalVcfPairInfo.pairId,
            name = "HighConfidence",
            bedpe = finalVcfPairInfo.svHighConfidenceFinalBedPe,
            chromLengths = chromLengths,
            listOfChroms = listOfChroms,
            diskSize = bedPeDiskSize
    }
    
    call tests.DescribeBedPe as allSomaticDescribeBedPe {
        input:
            pairId = finalVcfPairInfo.pairId,
            name = "AllSomatic",
            bedpe = finalVcfPairInfo.svFinalBedPe,
            chromLengths = chromLengths,
            listOfChroms = listOfChroms,
            diskSize = bedPeDiskSize
    }
    
    call tests.DraftSampleReport {
        input:
            listOfChroms = listOfChroms,
            karyotype = karyotype,
            chromLengths = chromLengths,
            docLink = docLink,
            
            pairId = finalVcfPairInfo.pairId,
            normal = finalVcfPairInfo.normal,
    
            cnvTable = DescribeBed.cnvTable,
            cnvGenesTable = DescribeBedGenes.cnvGenesTable,
    
            svGenesTable = DescribeBedPeGenes.svGenesTable,
            svTable = allSomaticDescribeBedPe.svTable,
            svHighConfidenceTable = highConfidenceDescribeBedPe.svTable,
    
            detailedGermVcfTable = SummarizeFinalGermVcf.detailedVcfTable,
            detailedVcfTable = SummarizeFinalVcf.detailedVcfTable,
            summaryVcfTable = SummarizeFinalVcf.summaryVcfTable,
            detailedLongOutputTable = SummarizeFinalVcf.detailedLongOutputTable,
            
            alleleCountsTxt = alleleCountsTxt,
            
            kouramiResult = kouramiResult,
            
            mantisStatusFinal = mantisStatusFinal,
            
            diff = diff,
            sig_input = sig_input,
            reconstructed = reconstructed,
            sigs = sigs,
            
            navTemplate=navTemplate,
            
            diskSize = 10
    }
    
    call tests.PrintReport {
        input:
            pairId = finalVcfPairInfo.pairId,
            md = DraftSampleReport.md,
            header = DraftSampleReport.header,
            navCustom = DraftSampleReport.navCustom,
            pandocTemplate = pandocTemplate,
            
            diskSize = 20
    }
    
    
    
    output {
        File cnvTable = DescribeBed.cnvTable
        File cnvGenesTable = DescribeBedGenes.cnvGenesTable
        
        File svGenesTable = DescribeBedPeGenes.svGenesTable
        File highConfidenceSvTable = highConfidenceDescribeBedPe.svTable
        File allSomaticSvTable = allSomaticDescribeBedPe.svTable
        
        File detailedGermVcfTable = SummarizeFinalGermVcf.detailedVcfTable
        File detailedVcfTable = SummarizeFinalVcf.detailedVcfTable
        File summaryVcfTable = SummarizeFinalVcf.summaryVcfTable
        
        File report = PrintReport.report

    }
}

