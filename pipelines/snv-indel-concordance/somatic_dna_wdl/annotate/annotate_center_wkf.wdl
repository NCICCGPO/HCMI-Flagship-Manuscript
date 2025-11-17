version 1.0

import "../merge_vcf/merge_vcf.wdl" as mergeVcf
import "annotate.wdl" as annotate
import "variantEffectPredictor.wdl" as variantEffectPredictor
import "../tasks/utils.wdl" as utils
import "../wdl_structs.wdl"

workflow Annotate {
    input {
        String tumorId
        String normalId
        String pairName
        File unannotatedVcf

        Boolean production = true

        String vepGenomeBuild
        Map[String, IndexedVcf] cosmicCodingChrom
        Map[String, IndexedVcf] cosmicNoncodingChrom
        Array[String]+ listOfChroms

        # Public
        Map[String, File] vepCacheChrom
        Map[String, File] annotationsChrom
        File plugins
        File dbNSFPReplacementLogic
        File MaxEntScan
        Map[String, IndexedTable] dbNSFP4Chrom
        Map[String, IndexedTable] dbscSNVChrom
        Map[String, IndexedVcf] gnomadGenomes
        Map[String, IndexedVcf] gnomadExomes
        IndexedReference vepFastaReference

        # Public
        File cancerResistanceMutations
        IndexedVcf deepIntronicsVcf
        IndexedVcf clinvarIntronicsVcf

        # post annotation
        File cosmicCensus

        String library

        IndexedReference referenceFa
        IndexedTable exonicRegionBedAnnotationsGz

    }

    # gather split VCF filenames to avoid glob that can fail on prem
    Array[String] suffixes = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29", "30", "31", "32", "33", "34", "35", "36", "37", "38", "39", "40", "41", "42", "43", "44", "45", "46", "47", "48", "49", "50"]
    scatter(chrom in listOfChroms) {
        Int maxSplits = 5
        String prefix = "~{pairName}.snv.indel.v7.~{chrom}"
        String additionalSuffix = ".vcf"
        Array[String] mnvVcfPaths = ["~{prefix}.mnv~{additionalSuffix}"]
        scatter (index in range(length(suffixes))) {
            String suffix = suffixes[index]
            String splitChromVcfPaths = "~{prefix}.~{suffix}~{additionalSuffix}"
        }
        Array[String] splitVcfPaths = flatten([mnvVcfPaths, splitChromVcfPaths])
        call utils.SplitVcfChrom {
            input:
                vcf = unannotatedVcf,
                prefix = prefix,
                diskSize = (ceil(size(unannotatedVcf, "GB")) * 3) + 10,
                chrom = chrom,
                maxRows = 1000,
                minSplits = 2,
                maxSplits = maxSplits,
                splitVcfPaths = splitVcfPaths
        }
        Array[File] splitVcfs = select_all(SplitVcfChrom.splitVcfs)
        scatter(vcf in splitVcfs) {
            call mergeVcf.CompressIndexVcf as unannotatedVcf {
                input:
                    vcf = vcf,
                    memoryGb = 8
            }
        
            Int vepDiskSize = ceil(size(vepCacheChrom[chrom], "GB") + size(gnomadGenomes[chrom].vcf, "GB") + size(gnomadExomes[chrom].vcf, "GB") + size(cosmicNoncodingChrom[chrom].vcf, "GB") + size(cosmicCodingChrom[chrom].vcf, "GB") + size(dbscSNVChrom[chrom].table, "GB") + size(dbNSFP4Chrom[chrom].table, "GB") + size(annotationsChrom[chrom], "GB") + size(deepIntronicsVcf.vcf, "GB") + size(clinvarIntronicsVcf.vcf, "GB") + (size(vcf, "GB") * 2)) + 20
            call variantEffectPredictor.vepPublicSvnIndel as productionVepSvnIndel {
                input:
                    pairName = pairName,
                    unannotatedVcf = unannotatedVcf.vcfCompressedIndexed,
                    vepCache = vepCacheChrom[chrom],
                    annotations = annotationsChrom[chrom],
                    plugins = plugins,
                    vepGenomeBuild = vepGenomeBuild,
                    vepFastaReference = vepFastaReference,
                    # Public
                    deepIntronicsVcf = deepIntronicsVcf,
                    clinvarIntronicsVcf = clinvarIntronicsVcf,
                    cosmicCoding = cosmicCodingChrom[chrom],
                    cosmicNoncoding = cosmicNoncodingChrom[chrom],
                    dbNSFPReplacementLogic = dbNSFPReplacementLogic,
                    MaxEntScan = MaxEntScan,
                    dbNSFP4 = dbNSFP4Chrom[chrom],
                    dbscSNV = dbscSNVChrom[chrom],
                    gnomadGenomes = gnomadGenomes[chrom],
                    gnomadExomes = gnomadExomes[chrom],
                    diskSize = vepDiskSize,
                    memoryGb = 8
            }
            call annotate.CustomAnnotateCenter {
                input:
                    pairName = pairName,
                    vcfAnnotatedVep = productionVepSvnIndel.vcfAnnotatedVep
            }
        }
    }
    Array[File] tempVcfs = flatten(CustomAnnotateCenter.vcfCsqRenamed)
    
    Int mergeDiskSize = (ceil(size(tempVcfs, "GB")) * 2) + 10
    call mergeVcf.Gatk4MergeSortVcf {
        input:
            tempVcfs = tempVcfs,
            sortedVcfPath = "~{pairName}.snv.indel.union.annotated.vcf",
            referenceFa = referenceFa,
            memoryGb = 8,
            diskSize = mergeDiskSize
    }
    
    if (library == "WGS") {
        call annotate.MainCenterVcf as wgsMainCenterVcf  {
        input:
            pairName = pairName,
            vcfAnnotated = Gatk4MergeSortVcf.sortedVcf.vcf
        }
    }
    if (library == "WgsExome"){
        call annotate.MainExonicCenterVcf as wgsExomeMainCenterVcf {
        input:
            pairName = pairName,
            exonicRegionBedAnnotationsGz = exonicRegionBedAnnotationsGz,
            vcfAnnotated = Gatk4MergeSortVcf.sortedVcf.vcf
        }
    }
    File runMainVcf = select_first([wgsMainCenterVcf.mainVcf, wgsExomeMainCenterVcf.mainVcf])
    File runMainIntersectionVcf = select_first([wgsMainCenterVcf.mainIntersectionVcf, wgsExomeMainCenterVcf.mainIntersectionVcf])

    call annotate.VcfToMultiCenterMaf {
        input:
            tumor = tumorId,
            normal = normalId,
            pairName = pairName,
            mainVcf = runMainVcf,
            intersectionMafPath = "~{pairName}.snv.indel.intersection.annotated.maf",
            unionMafPath = "~{pairName}.snv.indel.union.annotated.maf",
            library = library,
            vepGenomeBuild = vepGenomeBuild
    }

    output {
        File mainVcf = "~{runMainVcf}"
        File mainIntersectionVcf = "~{runMainIntersectionVcf}"
        File intersectionMaf = "~{VcfToMultiCenterMaf.intersectionMaf}"
        File unionMaf = "~{VcfToMultiCenterMaf.unionMaf}"
    }
}
