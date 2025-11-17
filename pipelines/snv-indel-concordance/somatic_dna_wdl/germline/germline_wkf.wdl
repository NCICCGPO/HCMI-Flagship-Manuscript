version 1.0

import "../wdl_structs.wdl"
import "germline.wdl" as germline
import "../merge_vcf/merge_vcf.wdl" as mergeVcf
import "../calling/calling.wdl" as calling

workflow Germline {
    # command

    input {
        Bam finalBam
        IndexedReference referenceFa
        Array[String]+ listOfChroms
        String normal
        String outputPrefix = "~{normal}"

        File excludeIntervalList
        Array[File] scatterIntervalsHcs

        IndexedVcf MillsAnd1000G
        IndexedVcf omni
        IndexedVcf hapmap
        IndexedReference referenceFa
        IndexedVcf onekG
        IndexedVcf dbsnp

        IndexedVcf whitelist
        IndexedVcf nygcAf
        IndexedVcf pgx
        IndexedTable rwgsPgxBed
        IndexedVcf deepIntronicsVcf
        IndexedVcf clinvarIntronicsVcf
        IndexedVcf chdWhitelistVcf

        Boolean highMem = false

        Int hcDiskSize = 20
    }
    scatter (i in range(length(scatterIntervalsHcs))) {
        call germline.haplotypeCallerGatk4 {
            input:
                referenceFa = referenceFa,
                finalBam = finalBam,
                index = i,
                outputPrefix = outputPrefix,
                excludeIntervalList = excludeIntervalList,
                scatterIntervalsHc = scatterIntervalsHcs[i],
                diskSize = hcDiskSize
        }
    }

    scatter (vcf in haplotypeCallerGatk4.haplotypecallerIntervalVcf) {
        File haplotypecallerIntervalVcf = vcf.vcf
    }
    Array[File] haplotypecallerIntervalVcfs = haplotypecallerIntervalVcf

    Int lowMergeSortDiskSize = 80

    if (highMem) {
        Int highMergeSortDiskSize = 200
    }

    Int mergeSortDiskSize = select_first([highMergeSortDiskSize, lowMergeSortDiskSize])

    call calling.Gatk4MergeSortVcf as haplotypecallerGatk4MergeSortVcf {
        input:
            sortedVcfPath = "~{outputPrefix}.haplotypecaller.g.vcf",
            tempChromVcfs = haplotypecallerIntervalVcfs,
            referenceFa = referenceFa,
            memoryGb = 16,
            diskSize = mergeSortDiskSize
    }

    call mergeVcf.CompressVcf as haplotypecallerCompressVcf {
        input:
            vcf = haplotypecallerGatk4MergeSortVcf.sortedVcf.vcf
    }

    call mergeVcf.IndexVcf as haplotypecallerIndexVcf {
        input:
            vcfCompressed = haplotypecallerCompressVcf.vcfCompressed
    }

    scatter (i in range(length(scatterIntervalsHcs))) {
        call germline.GentotypeGvcfsGatk4 {
            input:
                referenceFa = referenceFa,
                outputPrefix = outputPrefix,
                index = i,
                sortedVcf = haplotypecallerIndexVcf.vcfCompressedIndexed,
                scatterIntervalsHc=scatterIntervalsHcs[i],
                omni = omni,
                hapmap = hapmap,
                onekG = onekG,
                dbsnp = dbsnp
        }
    }

    scatter (vcf in GentotypeGvcfsGatk4.haplotypecallerFilteredGenoVcf) {
        File haplotypecallerFilteredGenoVcf = vcf.vcf
    }
    Array[File] haplotypecallerFilteredGenoVcfs = haplotypecallerFilteredGenoVcf

    call calling.Gatk4MergeSortVcf as genotypedFilteredMergeSortVcf {
        input:
            sortedVcfPath = "~{outputPrefix}.haplotypecaller.gatk.filtered.genotypedGVCFs.vcf",
            tempChromVcfs = haplotypecallerFilteredGenoVcfs,
            referenceFa = referenceFa,
            memoryGb = 16,
            diskSize = 200
    }

    call mergeVcf.CompressVcf as genotypedFilteredCompressVcf {
        input:
            vcf = genotypedFilteredMergeSortVcf.sortedVcf.vcf
    }

    call mergeVcf.IndexVcf as genotypedFilteredIndexVcf {
        input:
            vcfCompressed = genotypedFilteredCompressVcf.vcfCompressed
    }

    call germline.genotypeRefinementWorkflow {
        input:
            genotypedGatk4 = genotypedFilteredIndexVcf.vcfCompressedIndexed,
            outputPrefix = outputPrefix,
            referenceFa = referenceFa
    }

    Int lowHoDiskSize = 100
    if (highMem) {
        Int highHoDiskSize = 150
    }
    Int hoDiskSize = select_first([highHoDiskSize, lowHoDiskSize])

    call germline.filterHO as filterHO {
        input:
            sampleId = normal,
            outputPrefix = outputPrefix,
            nygcAf = nygcAf,
            haplotypecallerAfVcf = genotypeRefinementWorkflow.haplotypecallerAfVcf,
            pgx = pgx,
            rwgsPgxBed = rwgsPgxBed,
            whitelist = whitelist,
            chdWhitelistVcf = chdWhitelistVcf,
            deepIntronicsVcf = deepIntronicsVcf,
            clinvarIntronicsVcf = clinvarIntronicsVcf,
            diskSize = hoDiskSize
    }

    output {
        # equal to mergeHO.output_vcf
        IndexedVcf haplotypecallerVcf = genotypeRefinementWorkflow.haplotypecallerAfVcf
        # pass on to annotation
        IndexedVcf haplotypecallerFinalFiltered = filterHO.haplotypecallerFinalFiltered
    }
}
