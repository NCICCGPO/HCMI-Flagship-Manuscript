version 1.0

import "../merge_vcf/merge_vcf.wdl" as mergeVcf
import "annotate.wdl" as annotate
import "variantEffectPredictor.wdl" as variantEffectPredictor
import "../tasks/utils.wdl" as utils
import "../wdl_structs.wdl"

workflow GermlineAnnotate {
    input {
        String normal
        String sampleId = "~{normal}"
        IndexedVcf unannotatedVcf
        String haplotypecallerAnnotatedVcfPath = "~{normal}.haplotypecaller.gatk.annotated.vcf"
        Boolean production = true
        Boolean local = false

        Array[String]+ listOfChroms

        String vepGenomeBuild
        Map[String, IndexedVcf] cosmicCodingChrom
        Map[String, IndexedVcf] cosmicNoncodingChrom

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
        # split multi alleleic sites and remove HLA contig calls
        call mergeVcf.SplitMultiAllelicRegions {
                input:
                    pairName = sampleId,
                    vcfCompressedIndexed = unannotatedVcf,
                    listOfChroms = [chrom],
                    splitVcfPath = sub(basename(unannotatedVcf.vcf), ".vcf.gz$", ".split.vcf"),
                    referenceFa = referenceFa
            }
        String prefix = "~{sampleId}.split.haplotypecaller.gatk"
        String additionalSuffix = ".vcf"
        Array[String] mnvVcfPaths = ["~{prefix}.mnv~{additionalSuffix}"]
        scatter (index in range(length(suffixes))) {
            String suffix = suffixes[index]
            String splitChromVcfPaths = "~{prefix}.~{suffix}~{additionalSuffix}"
        }
        Array[String] splitVcfPaths = flatten([mnvVcfPaths, splitChromVcfPaths])
        call utils.SplitVcfChrom {
        input:
            vcf=SplitMultiAllelicRegions.sortedVcf,
            prefix=prefix,
            diskSize = (ceil(size(SplitMultiAllelicRegions.sortedVcf, "GB")) * 3) + 10,
            chrom = chrom,
            maxRows = 1000,
            minSplits = select_first([cloudMinSplits, hpcMinSplits]),
            maxSplits = select_first([cloudMaxSplits, hpcMaxSplits]),
            splitVcfPaths = splitVcfPaths
        }
        Array[File] splitVcfs = select_all(SplitVcfChrom.splitVcfs)
        scatter (vcf in splitVcfs) {
            Int vepDiskSize = ceil(size(vepCacheChrom[chrom], "GB") + size(gnomadGenomes[chrom].vcf, "GB") + size(gnomadExomes[chrom].vcf, "GB") + size(cosmicNoncodingChrom[chrom].vcf, "GB") + size(cosmicCodingChrom[chrom].vcf, "GB") + size(dbscSNVChrom[chrom].table, "GB") + size(dbNSFP4Chrom[chrom].table, "GB") + size(annotationsChrom[chrom], "GB") + size(deepIntronicsVcf.vcf, "GB") + size(clinvarIntronicsVcf.vcf, "GB") + (size(vcf, "GB") * 2)) + 20
            call mergeVcf.CompressIndexVcf as unannotatedVcf {
                input:
                    vcf = vcf,
                    memoryGb = 8
            }

            call variantEffectPredictor.vepPublicSvnIndel as productionVepSvnIndel {
                input:
                    pairName = sampleId,
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

            call annotate.RemoveSpanning {
                input:
                    sampleId = sampleId,
                    vcfAnnotatedVep = productionVepSvnIndel.vcfAnnotatedVep,
                    diskSize = ceil( size(productionVepSvnIndel.vcfAnnotatedVep, "GB") * 2) + 20
            }
            
            call annotate.CustomAnnotate {
                input:
                    pairName = sampleId,
                    cosmicCensus = cosmicCensus,
                    cancerResistanceMutations = cancerResistanceMutations,
                    vcfAnnotatedVep = RemoveSpanning.noSpanningVcf
            }
        }
    }
    
    Array[File] tempVcfs = flatten(CustomAnnotate.vcfCsqRenamed)
    
    Int mergeDiskSize = (ceil(size(tempVcfs, "GB")) * 3) + 10
    call mergeVcf.Gatk4MergeSortVcf {
        input:
            tempVcfs = tempVcfs,
            sortedVcfPath = "~{haplotypecallerAnnotatedVcfPath}",
            referenceFa = referenceFa,
            memoryGb = 32,
            diskSize = mergeDiskSize
    }

    output {

        File haplotypecallerAnnotatedVcf = Gatk4MergeSortVcf.sortedVcf.vcf
    }
}
