version 1.0

import "calling.wdl" as calling
import "../wdl_structs.wdl"

workflow Lancet {
    # command
    #   run Lancet caller
    input {
        Boolean local
        String library
        String tumorId
        String normalId
        Array[String]+ listOfChroms
        Array[File] splitBedsWgs
        Array[File] splitBedsWgsFewNodes
        Map[String, File] chromBeds
        String pairName
        IndexedReference referenceFa
        Bam normalFinalBam
        Bam tumorFinalBam
        # resources
        Int threads = 8
        Int memoryGb = 40
        File lancetJsonLog
        
        # needed for work on cloud without localizing
        File? serviceAccountKey
        String? gcpProject

    }

    if (library == 'Exome') {
        scatter(chrom in listOfChroms) {
            File chromBedsExome = chromBeds[chrom]
        }
    }
    if (local) {
        Array[File] splitBedsWgsRun = splitBedsWgsFewNodes
    }
    Array[File] splitBeds = select_first([splitBedsWgsRun, splitBedsWgs, chromBedsExome])
    
    scatter(splitBed in splitBeds) {
        call calling.SliceBam as sliceBamTumor {
            input:
                chrom = 'split',
                sampleId = tumorId,
                chromBed = splitBed,
                finalBam = tumorFinalBam,
                threads = 4,
                memoryGb = 4,
                diskSize = 10,
                serviceAccountKey = serviceAccountKey,
                gcpProject = gcpProject
        }
        
        call calling.SliceBam as sliceBamNormal {
            input:
                chrom = 'split',
                sampleId = normalId,
                chromBed = splitBed,
                finalBam = normalFinalBam,
                threads = 4,
                memoryGb = 4,
                diskSize = 10,
                serviceAccountKey = serviceAccountKey,
                gcpProject = gcpProject
        }
        
        Int diskSize = ceil( size(sliceBamNormal.lancetBamSlice.bam, "GB") + size(sliceBamTumor.lancetBamSlice.bam, "GB")) + 10
                
        call calling.LancetWGSRegional {
            input:
                chrom = 'split',
                lancetChromVcfPath = "~{pairName}.split.lancet.vcf",
                chromBed = splitBed,
                referenceFa = referenceFa,
                normalFinalBam = sliceBamNormal.lancetBamSlice,
                tumorFinalBam = sliceBamTumor.lancetBamSlice,
                pairName = pairName,
                threads = threads,
                memoryGb = memoryGb,
                diskSize = diskSize
        }
    }
    
    call calling.Gatk4MergeSortVcf {
        input:
            sortedVcfPath = "~{pairName}.lancet.sorted.vcf",
            tempChromVcfs = LancetWGSRegional.lancetChromVcf,
            referenceFa = referenceFa,
            memoryGb = 8,
            diskSize = 10
    }

    call calling.AddVcfCommand as lancetAddVcfCommand {
        input:
            inVcf = Gatk4MergeSortVcf.sortedVcf.vcf,
            jsonLog = lancetJsonLog,
            memoryGb = 4,
            diskSize = 10
    }

    call calling.ReorderVcfColumns as lancetReorderVcfColumns {
        input:
            tumor = tumorId,
            normal = normalId,
            rawVcf = lancetAddVcfCommand.outVcf,
            orderedVcfPath = "~{pairName}.lancet.vcf",
            memoryGb = 4,
            diskSize = 10
    }

    output {
        File lancet = lancetReorderVcfColumns.orderedVcf
    }
}
