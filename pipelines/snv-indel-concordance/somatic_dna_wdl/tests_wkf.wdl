version 1.0

import "test/tests_wkf.wdl" as tests
import "wdl_structs.wdl"

# ================== COPYRIGHT ================================================
# New York Genome Center
# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright (2021) by the New York
# Genome Center. All rights are reserved. This software is supplied without
# any warranty or guaranteed support whatsoever. The New York Genome Center
# cannot be responsible for its use, misuse, or functionality.
#
#    Jennifer M Shelton (jshelton@nygenome.org)
#    Nico Robine (nrobine@nygenome.org)
#    Minita Shah (mshah@nygenome.org)
#    Timothy Chu (tchu@nygenome.org)
#    Will Hooper (whooper@nygenome.org)
#
# ================== /COPYRIGHT ===============================================

workflow PipelineTests {
    input {
        File chromLengths
        File cosmicCensus
        Array[String] listOfChroms
        # BAMs
        Array[File] newflagStat
        Array[File] oldflagStat
        # SNV/INDEL calls
        Array[PairVcfInfo] oldpairVcfInfos
        Array[PairVcfInfo] newpairVcfInfos
        # CNV/SV calls
        Array[File] oldcnvAnnotatedFinalBed
        Array[File] newcnvAnnotatedFinalBed
        Array[File] oldsvFinalBedPe
        Array[File] newsvFinalBedPe
        # MSI
        Array[File] oldmantisStatusFinal
        Array[File] newmantisStatusFinal
        # baf
        Array[File] oldalleleCountsTxt
        Array[File] newalleleCountsTxt
        # Hla
        Array[File] oldkouramiResult
        Array[File] newkouramiResult
        IndexedReference referenceFa
    }
    
    scatter (pairVcfInfo in oldpairVcfInfos) {
        String pairIds = pairVcfInfo.pairId
    }
    
    call tests.Tests {
        input:
            chromLengths = chromLengths,
            cosmicCensus = cosmicCensus,
            listOfChroms = listOfChroms,
            oldBafs = oldalleleCountsTxt,
            newBafs = newalleleCountsTxt,
            
            newflagStat = newflagStat,
            oldflagStat = oldflagStat,
            
            oldmantisStatusFinal = oldmantisStatusFinal,
            newmantisStatusFinal = newmantisStatusFinal,
            
            oldkouramiResult = oldkouramiResult,
            newkouramiResult = newkouramiResult,
            
            pairIds = pairIds,
            
            oldpairVcfInfos = oldpairVcfInfos,
            newpairVcfInfos = newpairVcfInfos,
            oldcnvAnnotatedFinalBed = oldcnvAnnotatedFinalBed,
            newcnvAnnotatedFinalBed = newcnvAnnotatedFinalBed,
            oldsvFinalBedPe = oldsvFinalBedPe,
            newsvFinalBedPe = newsvFinalBedPe,
            referenceFa = referenceFa
    }

    output {
        File finalCnvGenesTable = Tests.finalCnvGenesTable
        File finalSvGenesTable = Tests.finalSvGenesTable
        File flagstatOutputTable = Tests.flagstatOutputTable
        File msiOutputTable = Tests.msiOutputTable
        File hlaOutputTable = Tests.hlaOutputTable
        File hlaOutputConcordance = Tests.hlaOutputConcordance
        File finalBedTable = Tests.finalBedTable
        File finalBedPeTable = Tests.finalBedPeTable
        File finalVcfTable = Tests.finalVcfTable
        Array[File] finalBedpeJoined = Tests.finalBedpeJoined
        Array[File] finalBedJoined = Tests.finalBedJoined
        Array[File] detailedOutputVcfTable = Tests.detailedOutputVcfTable
        File md = Tests.md
        File header = Tests.header
    }
    
    meta {
        usage : "Run on output from v7 all lists of output must be in the same order (based on the sampleId/pairId scatter used to create the output)"
        notes : "Append new or old to the begining of the WDL name found in the cromwell log or the custom outputInfo.json. For the moment sampleId order is assumed from filenames (needs updating)"
    }
}