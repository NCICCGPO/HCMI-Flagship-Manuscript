version 1.0

import "alignment_analysis.wdl" as alignmentAnalysis
import "../wdl_structs.wdl"

workflow Kourami {
    input {
        String sampleId
        String analysisId
        # mergedHlaPanel
        BwaReference kouramiReference
        IndexedReference referenceFa
        File kouramiFastaGem1Index
        Bam finalBam
        Int diskSize = ceil( size(finalBam.bam, "GB")) + 20
        Int refDiskSize = ceil( size(referenceFa.fasta, "GB")) + 1
        
        # needed for work on cloud without localizing
        File? serviceAccountKey
        String? gcpProject
    }
    
    call alignmentAnalysis.GetChr6Contigs {
        input:
            referenceFa = referenceFa,
            diskSize = refDiskSize
    }
    call alignmentAnalysis.GemSelect {
        input:
            sampleId = sampleId,
            finalBam = finalBam,
            chr6Contigs = GetChr6Contigs.chr6Contigs,
            kouramiFastaGem1Index = kouramiFastaGem1Index,
            diskSize = diskSize
    }
    
    call alignmentAnalysis.LookUpMates {
        input:
            sampleId = sampleId,
            r2File = GemSelect.r2File,
            r2MappedFastq = GemSelect.r2MappedFastq,
            r1File = GemSelect.r1File,
            r1MappedFastq = GemSelect.r1MappedFastq
    }
    
    if (defined(serviceAccountKey)) {
        call alignmentAnalysis.GetMates {
            input:
                sampleId = sampleId,
                finalBam = finalBam,
                r1UnmappedFile = LookUpMates.r1UnmappedFile,
                r2UnmappedFile = LookUpMates.r2UnmappedFile,
                diskSize = 20,
                serviceAccountKey = serviceAccountKey,
                gcpProject = gcpProject
        }
    }
    
    if (!defined(serviceAccountKey)) {
        call alignmentAnalysis.GetMatesLocalize {
            input:
                sampleId = sampleId,
                finalBam = finalBam,
                r1UnmappedFile = LookUpMates.r1UnmappedFile,
                r2UnmappedFile = LookUpMates.r2UnmappedFile,
                diskSize = diskSize
        }
    }
    
    File r1UnmappedFastq = select_first([GetMates.r1UnmappedFastq, GetMatesLocalize.r1UnmappedFastq])
    File r2UnmappedFastq = select_first([GetMates.r2UnmappedFastq, GetMatesLocalize.r2UnmappedFastq])
    String sortedR1FastqPath = "~{analysisId}.R1_sorted.fastq"
    call alignmentAnalysis.SortFastqs as r1SortFastqs {
        input:
            fastqPairId = "R1",
            sampleId = sampleId,
            sortedFastqPath = sortedR1FastqPath,
            chr6MappedFastq = GemSelect.r1MappedFastq,
            chr6MappedMatesFastq = r1UnmappedFastq
    }
    String sortedR2FastqPath = "~{analysisId}.R2_sorted.fastq"
    call alignmentAnalysis.SortFastqs as r2SortFastqs {
        input:
            fastqPairId = "R2",
            sampleId = sampleId,
            sortedFastqPath = sortedR2FastqPath,
            chr6MappedFastq = GemSelect.r2MappedFastq,
            chr6MappedMatesFastq = r2UnmappedFastq
    }
    
    call alignmentAnalysis.AlignToPanel {
        input:
            sampleId = sampleId,
            kouramiReference = kouramiReference,
            r1SortedFastq = r1SortFastqs.sortedFastq,
            r2SortedFastq = r2SortFastqs.sortedFastq
    }
    String resultPrefix = "~{analysisId}.kourami"
    call alignmentAnalysis.Kourami {
        input:
            sampleId = sampleId,
            resultPrefix = resultPrefix,
            kouramiBam = AlignToPanel.kouramiBam
    }
    
    output {
        File result = Kourami.result
        File r1HlaFastq = r1SortFastqs.sortedFastq
        File r2HlaFastq = r2SortFastqs.sortedFastq
    }
}
