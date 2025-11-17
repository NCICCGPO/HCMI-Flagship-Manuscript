version 1.0

import "tests.wdl" as tests
import "../tasks/utils.wdl" as utils
import "../wdl_structs.wdl"


workflow Tests {
    input {
        String name = "comparison"
        File chromLengths
        File cosmicCensus
        Array[String] listOfChroms
        # BAM
        Array[File] newflagStat
        Array[File] oldflagStat
        
        Array[String] pairIds
        
        Array[PairVcfInfo] oldpairVcfInfos
        Array[PairVcfInfo] newpairVcfInfos
        Array[File] oldcnvAnnotatedFinalBed
        Array[File] newcnvAnnotatedFinalBed
        Array[File] oldsvFinalBedPe
        Array[File] newsvFinalBedPe
        Array[File] oldBafs
        Array[File] newBafs
        # MSI
        Array[File] oldmantisStatusFinal
        Array[File] newmantisStatusFinal
        # Hla
        Array[File] oldkouramiResult
        Array[File] newkouramiResult
        IndexedReference referenceFa
    }
    
    scatter (flagStat in zip(newflagStat, oldflagStat)) {
        String sampleIdFlagstat = basename(flagStat.left, ".FlagStat.txt") 
        call tests.SummarizeFlagstat {
            input:
                sampleId = sampleIdFlagstat,
                newflagStat = flagStat.left,
                oldflagStat = flagStat.right
        }
    }
    
    scatter (mantisStatusFinal in zip(newmantisStatusFinal, oldmantisStatusFinal)) {

        String sampleIdMantis = basename(mantisStatusFinal.left, ".mantis.WGS-targeted.status.final.tsv") 
        call tests.SummarizeMsi {
            input:
                sampleId = sampleIdMantis,
                newmantisStatusFinal = mantisStatusFinal.left,
                oldmantisStatusFinal = mantisStatusFinal.right
        }
    }
    
    scatter (kouramiResult in zip(newkouramiResult, oldkouramiResult)) {
        String sampleIdKourami = basename(kouramiResult.left, ".kourami.result") 
        call tests.SummarizeHla {
            input:
                sampleId = sampleIdKourami,
                newkouramiResult = kouramiResult.left,
                oldkouramiResult = kouramiResult.right
        }
    }
    
    scatter (pairId in pairIds) {
        call utils.GetIndex {
            input:
                sampleIds = pairIds,
                sampleId = pairId
        }
        
        call tests.CompareBedPe as finalCompareBedPe {
            input:
                pairId = pairId,
                oldBedpe = oldsvFinalBedPe[GetIndex.index],
                newBedpe = newsvFinalBedPe[GetIndex.index]
        }
        
        call tests.CompareSvGenes as finalCompareSvGenes {
            input:
                pairId = pairId,
                name = name,
                concordanceBedPe = finalCompareBedPe.outFileBedpe
        }
        
        call tests.CompareBed as finalCompareBed {
            input:
                pairId = pairId,
                oldBed = oldcnvAnnotatedFinalBed[GetIndex.index],
                newBed = newcnvAnnotatedFinalBed[GetIndex.index]
        }
        
        call tests.CompareCnvGenes as finalCompareCnvGenes {
            input:
                pairId = pairId,
                name = name,
                concordanceBed = finalCompareBed.outFileBed
        }
        
        call tests.SomPy {
            input:
                referenceFa = referenceFa,
                pairId = pairId,
                oldVcf = oldpairVcfInfos[GetIndex.index].mainVcf,
                newVcf = newpairVcfInfos[GetIndex.index].mainVcf
        }
        
        call tests.SummarizeVcf as finalSummarizeVcf {
            input:
                pairId = pairId,
                name = "final.vcf",
                oldOnlyVcf = SomPy.oldOnlyVcf,
                newOnlyVcf = SomPy.newOnlyVcf,
                concordantVcf = SomPy.concordantVcf
        }
    }
    
    call tests.ConcateTables as finalCnvGenesConcateTables {
        input:
            tables = finalCompareCnvGenes.cnvGenesTable,
            outputTablePath = "~{name}.finalCnvGenes.csv"
    }
    
    call tests.ConcateTables as finalSvGenesConcateTables {
        input:
            tables = finalCompareSvGenes.svGenesTable,
            outputTablePath = "~{name}.finalSvGenes.csv"
    }
    
    call tests.ConcateTables as flagstatConcateTables {
        input:
            tables = SummarizeFlagstat.outputTable,
            outputTablePath = "~{name}.flagstat.csv"
    }
    
    call tests.ConcateTables as msiConcateTables {
        input:
            tables = SummarizeMsi.outputTable,
            outputTablePath = "~{name}.msi.csv"
    }
    
    call tests.ConcateTables as hlaConcateTables {
        input:
            tables = SummarizeHla.outputTable,
            outputTablePath = "~{name}.hla.csv"
    }
    
    call tests.ConcateTables as hlaConcordanceConcateTables {
        input:
            tables = SummarizeHla.outputConcordance,
            outputTablePath = "~{name}.hla.concordance.csv"
    }
    
    call tests.ConcateTables as finalBedConcateTables {
        input:
            tables = finalCompareBed.outFileSummary,
            outputTablePath = "~{name}.finalBed.tsv"
    }
    
    call tests.ConcateTables as finalBedPeConcateTables {
        input:
            tables = finalCompareBedPe.outFileSummary,
            outputTablePath = "~{name}.finalBedPe.tsv"
    }
    
    call tests.ConcateTables as finalVcfConcateTables {
        input:
            tables = finalSummarizeVcf.summaryOutputTable,
            outputTablePath = "~{name}.finalVcf.csv"
    }
    
    call tests.DraftComparison as finalDraftComparison {
        input:
            name = name,
            chromLengths = chromLengths,
            cosmicCensus = cosmicCensus,
            listOfChroms = listOfChroms,
            oldBafs = oldBafs,
            newBafs = newBafs,
            cnvBeds = finalCompareBed.outFileBed,
            svBedPes = finalCompareBedPe.outFileBedpe,
            vcfDetails = finalSummarizeVcf.detailedOutputTable,
            bedCounts = finalBedConcateTables.outputTable,
            bedPeCounts = finalBedPeConcateTables.outputTable,
            bedPeGenes = finalSvGenesConcateTables.outputTable,
            bedGenes = finalCnvGenesConcateTables.outputTable,
            vcfCounts = finalVcfConcateTables.outputTable,
            flagstatSummary = flagstatConcateTables.outputTable,
            msiSummary = msiConcateTables.outputTable
    }
    
    output {
        File finalCnvGenesTable = finalCnvGenesConcateTables.outputTable
        File finalSvGenesTable = finalSvGenesConcateTables.outputTable
        File flagstatOutputTable = flagstatConcateTables.outputTable
        File msiOutputTable = msiConcateTables.outputTable
        File hlaOutputTable = hlaConcateTables.outputTable
        File hlaOutputConcordance = hlaConcordanceConcateTables.outputTable
        File finalBedTable = finalBedConcateTables.outputTable
        File finalBedPeTable = finalBedPeConcateTables.outputTable
        File finalVcfTable = finalVcfConcateTables.outputTable
        Array[File] finalBedpeJoined = finalCompareBedPe.outFileBedpe
        Array[File] finalBedJoined = finalCompareBed.outFileBed
        Array[File] detailedOutputVcfTable = finalSummarizeVcf.detailedOutputTable
        File md = finalDraftComparison.md
        File header = finalDraftComparison.header

    }
}

