version 1.0

import "calling.wdl" as calling
import "../merge_vcf/merge_vcf.wdl" as mergeVcf
import "../wdl_structs.wdl"

workflow Gridss {
    # command
    #   run Gridss caller
    input {
        String tumor
        String normal
        String pairName

        BwaReference bwaReference
        Array[File] gridssAdditionalReference
        Array[String]+ listOfChroms

        Bam normalFinalBam
        Bam tumorFinalBam

        Int assembleChunks = 5

        String bsGenome
        File ponTarGz

        Int threads = 2

        File? gridssUnfilteredVcfChromsPreexist
        Boolean highMem = false
        Int preMemoryGb = 60
        Int filterMemoryGb = 32
    }

    Int tumorDiskSize = ceil(size(tumorFinalBam.bam, "GB") * 3) + 100
    Int normalDiskSize = ceil(size(normalFinalBam.bam, "GB") * 3) + 20
    if (!defined(gridssUnfilteredVcfChromsPreexist)) {

        call calling.GridssPreprocess as tumorGridssPreprocess {
            input:
                threads = threads,
                memoryGb = preMemoryGb,
                bwaReference = bwaReference,
                gridssAdditionalReference = gridssAdditionalReference,
                finalBam = tumorFinalBam,
                diskSize = tumorDiskSize
        }
    
        call calling.GridssPreprocess as normalGridssPreprocess {
            input:
                threads = threads,
                memoryGb = preMemoryGb,
                bwaReference= bwaReference,
                gridssAdditionalReference = gridssAdditionalReference,
                finalBam = normalFinalBam,
                diskSize = normalDiskSize
        }
    
        Int lowAssembleMemoryGb = 60
        Int lowAssembleDiskSize = ceil( size(tumorFinalBam.bam, "GB") * 1.4 ) + ceil( size(normalFinalBam.bam, "GB")  * 1.4)
    
        if (highMem) {
            Int highAssembleMemoryGb = 100
            Int highAssembleDiskSize = ceil( size(tumorFinalBam.bam, "GB") * 1.4 ) + ceil( size(normalFinalBam.bam, "GB")  * 1.4) + 20
        }
    
        Int assembleMemoryGb = select_first([highAssembleMemoryGb, lowAssembleMemoryGb])
        Int assembleDiskSize = select_first([highAssembleDiskSize, lowAssembleDiskSize])
        scatter(i in range(assembleChunks)) {
            call calling.GridssAssembleChunk {
                input:
                    threads = threads,
                    memoryGb = assembleMemoryGb,
                    pairName = pairName,
                    bwaReference = bwaReference,
                    gridssAdditionalReference = gridssAdditionalReference,
                    jobIndex = i,
                    assembleChunks = assembleChunks,
                    tumorFinalBam = tumorFinalBam,
                    normalFinalBam = normalFinalBam,
                    normalSvBam = normalGridssPreprocess.svBam,
                    normalCigarMetrics = normalGridssPreprocess.cigarMetrics,
                    normalIdsvMetrics = normalGridssPreprocess.idsvMetrics,
                    normalTagMetrics = normalGridssPreprocess.tagMetrics,
                    normalMapqMetrics = normalGridssPreprocess.mapqMetrics,
                    normalInsertSizeMetrics = normalGridssPreprocess.insertSizeMetrics,
                    tumorSvBam = tumorGridssPreprocess.svBam,
                    tumorCigarMetrics = tumorGridssPreprocess.cigarMetrics,
                    tumorIdsvMetrics = tumorGridssPreprocess.idsvMetrics,
                    tumorTagMetrics = tumorGridssPreprocess.tagMetrics,
                    tumorMapqMetrics = tumorGridssPreprocess.mapqMetrics,
                    tumorInsertSizeMetrics = tumorGridssPreprocess.insertSizeMetrics,
                    diskSize = assembleDiskSize
            }
        }

        call calling.GridssAssemble {
            input:
                threads = threads,
                memoryGb = assembleMemoryGb,
                pairName = pairName,
                bwaReference = bwaReference,
                gridssAdditionalReference = gridssAdditionalReference,
                tumorFinalBam = tumorFinalBam,
                normalFinalBam = normalFinalBam,
                downsampled = GridssAssembleChunk.downsampled,
                excluded = GridssAssembleChunk.excluded,
                subsetCalled = GridssAssembleChunk.subsetCalled,
                normalSvBam = normalGridssPreprocess.svBam,
                normalCigarMetrics = normalGridssPreprocess.cigarMetrics,
                normalIdsvMetrics = normalGridssPreprocess.idsvMetrics,
                normalTagMetrics = normalGridssPreprocess.tagMetrics,
                normalMapqMetrics = normalGridssPreprocess.mapqMetrics,
                normalInsertSizeMetrics = normalGridssPreprocess.insertSizeMetrics,
                tumorSvBam = tumorGridssPreprocess.svBam,
                tumorCigarMetrics = tumorGridssPreprocess.cigarMetrics,
                tumorIdsvMetrics = tumorGridssPreprocess.idsvMetrics,
                tumorTagMetrics = tumorGridssPreprocess.tagMetrics,
                tumorMapqMetrics = tumorGridssPreprocess.mapqMetrics,
                tumorInsertSizeMetrics = tumorGridssPreprocess.insertSizeMetrics,
                diskSize = assembleDiskSize
    
        }
    
        call calling.GridssCalling {
            input:
                threads = threads,
                memoryGb = assembleMemoryGb,
                pairName = pairName,
                bwaReference = bwaReference,
                gridssAdditionalReference = gridssAdditionalReference,
                tumorFinalBam = tumorFinalBam,
                normalFinalBam = normalFinalBam,
                downsampled = GridssAssembleChunk.downsampled,
                excluded = GridssAssembleChunk.excluded,
                subsetCalled = GridssAssembleChunk.subsetCalled,
                gridssassemblyBam = GridssAssemble.gridssassemblyBam,
                gridssassemblySvBam = GridssAssemble.gridssassemblySvBam,
                normalSvBam = normalGridssPreprocess.svBam,
                normalCigarMetrics = normalGridssPreprocess.cigarMetrics,
                normalIdsvMetrics = normalGridssPreprocess.idsvMetrics,
                normalTagMetrics = normalGridssPreprocess.tagMetrics,
                normalMapqMetrics = normalGridssPreprocess.mapqMetrics,
                normalInsertSizeMetrics = normalGridssPreprocess.insertSizeMetrics,
                tumorSvBam = tumorGridssPreprocess.svBam,
                tumorCigarMetrics = tumorGridssPreprocess.cigarMetrics,
                tumorIdsvMetrics = tumorGridssPreprocess.idsvMetrics,
                tumorTagMetrics = tumorGridssPreprocess.tagMetrics,
                tumorMapqMetrics = tumorGridssPreprocess.mapqMetrics,
                tumorInsertSizeMetrics = tumorGridssPreprocess.insertSizeMetrics,
                diskSize = assembleDiskSize + 20
    
        }
        Int lowNonChromsDiskSize = ceil( size(GridssCalling.gridssUnfilteredVcf, "GB")) + 30
        if (highMem) {
            Int highNonChromsDiskSize = ceil( size(GridssCalling.gridssUnfilteredVcf, "GB")) + 60
        }
        Int filterNonChromsDiskSize = select_first([highNonChromsDiskSize, lowNonChromsDiskSize])
        call calling.FilterNonChroms {
            input:
                diskSize = filterNonChromsDiskSize,
                memoryGb = 1,
                pairName = pairName,
                gridssUnfilteredVcf = GridssCalling.gridssUnfilteredVcf,
                listOfChroms = listOfChroms
        }
    }
    
    File gridssUnfilteredVcfChromsRun = select_first([FilterNonChroms.gridssUnfilteredVcfChroms, gridssUnfilteredVcfChromsPreexist])

    Int lowFilterDiskSize = ceil( size(gridssUnfilteredVcfChromsRun, "GB")) + 30
    if (highMem) {
        Int highFilterDiskSize = ceil( size(gridssUnfilteredVcfChromsRun, "GB")) + 60
    }
    Int filterDiskSize = select_first([highFilterDiskSize, lowFilterDiskSize])

    call calling.GridssFilter {
        input:
            threads = threads,
            memoryGb = filterMemoryGb,
            pairName = pairName,
            bsGenome = bsGenome,
            ponTarGz = ponTarGz,
            gridssUnfilteredVcf = gridssUnfilteredVcfChromsRun,
            diskSize = filterDiskSize
    }

    call calling.ReorderVcfColumns as filteredReorderVcfColumns {
        input:
            tumor = tumor,
            normal = normal,
            rawVcf = GridssFilter.gridssVcf.vcf,
            orderedVcfPath = "~{pairName}.sv.gridss.vcf",
            memoryGb = 4,
            diskSize = 10
    }

    call mergeVcf.CompressVcf {
        input:
            vcf = filteredReorderVcfColumns.orderedVcf,
            memoryGb = 4
    }

    call mergeVcf.IndexVcf {
        input:
            vcfCompressed = CompressVcf.vcfCompressed
    }

    output {
        IndexedVcf gridssVcf = IndexVcf.vcfCompressedIndexed
    }
}
