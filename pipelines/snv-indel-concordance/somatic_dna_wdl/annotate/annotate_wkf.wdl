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
        Boolean local = false

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

        File ensemblEntrez
        String library

        IndexedReference referenceFa

    }
    if (!local) {
        Int cloudMaxSplits = 5
        Int cloudMinSplits = 2
    }
    if (local) {
        Int hpcMaxSplits = 1
        Int hpcMinSplits = 1
    }
    
    # gather split VCF filenames to avoid glob that can fail on prem
    Array[String] suffixes = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29", "30", "31", "32", "33", "34", "35", "36", "37", "38", "39", "40", "41", "42", "43", "44", "45", "46", "47", "48", "49", "50"]
    scatter(chrom in listOfChroms) {
#         Int maxSplits = 5
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
                minSplits = select_first([cloudMinSplits, hpcMinSplits]),
                maxSplits = select_first([cloudMaxSplits, hpcMaxSplits]),
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
                    memoryGb = 7
            }
            
            call annotate.CustomAnnotate {
                input:
                    pairName = pairName,
                    cosmicCensus = cosmicCensus,
                    cancerResistanceMutations = cancerResistanceMutations,
                    vcfAnnotatedVep = productionVepSvnIndel.vcfAnnotatedVep
            }
        }
    }
    Array[File] tempVcfs = flatten(CustomAnnotate.vcfCsqRenamed)
    
    Int mergeDiskSize = (ceil(size(tempVcfs, "GB")) * 2) + 10
    call mergeVcf.Gatk4MergeSortVcf {
        input:
            tempVcfs = tempVcfs,
            sortedVcfPath = "~{pairName}.snv.indel.supplemental.v7.annotated.vcf",
            referenceFa = referenceFa,
            memoryGb = 8,
            diskSize = mergeDiskSize
    }

    call annotate.MainVcf {
        input:
            pairName = pairName,
            vcfAnnotated = Gatk4MergeSortVcf.sortedVcf.vcf
    }

    call annotate.TableVcf {
        input:
            tumor = tumorId,
            normal = normalId,
            pairName = pairName,
            mainVcf = MainVcf.mainVcf
    }

    call annotate.VcfToMaf {
        input:
            tumor = tumorId,
            normal = normalId,
            pairName = pairName,
            mainVcf = MainVcf.mainVcf,
            library = library,
            vepGenomeBuild = vepGenomeBuild,
            ensemblEntrez = ensemblEntrez
    }

    output {
        PairVcfInfo pairVcfInfo  = object {
            pairId : "~{pairName}",
            tumor : "~{tumorId}",
            normal : "~{normalId}",
            mainVcf : "~{MainVcf.mainVcf}",
            supplementalVcf : "~{Gatk4MergeSortVcf.sortedVcf.vcf}",
            vcfAnnotatedTxt : "~{TableVcf.vcfAnnotatedTxt}",
            maf : "~{VcfToMaf.maf}",
        }
    }
}
