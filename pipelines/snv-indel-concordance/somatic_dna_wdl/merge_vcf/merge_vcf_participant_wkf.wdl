version 1.0

import "merge_aliquots_wkf.wdl" as mergeAliquots
import "merge_vcf.wdl" as mergeVcf
import "../annotate/variantEffectPredictor.wdl" as variantEffectPredictor
import "../tasks/utils.wdl" as utils
import "../wdl_structs.wdl"

workflow MergeParticipantVcf {
    # command
    #   Call variants in BAMs
    #   merge and filter raw VCFs
    #   annotate
    input {
        Boolean local = false
        Array[File] centerVcfs
        Array[String] tumorBarcodeAliquots
        Array[String] nonNormalKinds
        String normalBarcodeAliquot
        Array[Bam] tumorFinalBams
        Bam normalFinalBam
        String participantId
        File cosmicCensus
        
        IndexedReference referenceFa
        Array[String]+ listOfChroms
        
        # optional annotation
        Boolean reAnnotate = false
        String vepGenomeBuild
        Map[String, IndexedVcf] cosmicCodingChrom
        Map[String, IndexedVcf] cosmicNoncodingChrom
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
        File cancerResistanceMutations
        IndexedVcf deepIntronicsVcf
        IndexedVcf clinvarIntronicsVcf
        # merge callers
        File? serviceAccountKey
        String? gcpProject
    }
    
    scatter (i in range(length(centerVcfs))) {
        if (i > 0) {
            # keep normal in first VCF
            call mergeVcf.SelectVcfSamples {
                input:
                    retainedSamples = [tumorBarcodeAliquots[i]],
                    centerVcf = centerVcfs[i],
                    referenceFa = referenceFa
            }
        }
        File runCenterVcf = select_first([SelectVcfSamples.nonNormalVcf, centerVcfs[i]])
        if (reAnnotate) {
            call mergeVcf.WipeAnnnotation {
                input:
                    inputVcf = runCenterVcf,
                    outputVcfPath = "~{tumorBarcodeAliquots[i]}--~{normalBarcodeAliquot}.csq_wiped.vcf"
            }
            
            Array[String] suffixes = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29", "30", "31", "32", "33", "34", "35", "36", "37", "38", "39", "40", "41", "42", "43", "44", "45", "46", "47", "48", "49", "50"]
            scatter(chrom in listOfChroms) {
                Int maxSplits = 5
                String pairName = "~{tumorBarcodeAliquots[i]}--~{normalBarcodeAliquot}"
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
                        vcf = WipeAnnnotation.outputVcf,
                        prefix = prefix,
                        diskSize = (ceil(size(WipeAnnnotation.outputVcf, "GB")) * 3) + 10,
                        chrom = chrom,
                        maxRows = 1000,
                        minSplits = 2,
                        maxSplits = maxSplits,
                        splitVcfPaths = splitVcfPaths
                }
                Array[File] splitVcfs = select_all(SplitVcfChrom.splitVcfs)
                scatter(vcf in splitVcfs) {
                    call mergeVcf.CompressIndexVcf as unannotatedCompressIndexVcf {
                        input:
                            vcf = vcf,
                            memoryGb = 8
                    }
                
                    Int vepDiskSize = ceil(size(vepCacheChrom[chrom], "GB") + size(gnomadGenomes[chrom].vcf, "GB") + size(gnomadExomes[chrom].vcf, "GB") + size(cosmicNoncodingChrom[chrom].vcf, "GB") + size(cosmicCodingChrom[chrom].vcf, "GB") + size(dbscSNVChrom[chrom].table, "GB") + size(dbNSFP4Chrom[chrom].table, "GB") + size(annotationsChrom[chrom], "GB") + size(deepIntronicsVcf.vcf, "GB") + size(clinvarIntronicsVcf.vcf, "GB") + (size(vcf, "GB") * 2)) + 20
                    call variantEffectPredictor.vepPublicSvnIndel as productionVepSvnIndel {
                        input:
                            pairName = pairName,
                            unannotatedVcf = unannotatedCompressIndexVcf.vcfCompressedIndexed,
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
                }
            }
            Array[File] tempVcfs = flatten(productionVepSvnIndel.vcfAnnotatedVep)
            Int mergeDiskSize = (ceil(size(tempVcfs, "GB")) * 2) + 10
            call mergeVcf.Gatk4MergeSortVcf {
                input:
                    tempVcfs = tempVcfs,
                    sortedVcfPath = "~{pairName}.snv.indel.union.annotated.vcf",
                    referenceFa = referenceFa,
                    memoryGb = 8,
                    diskSize = mergeDiskSize
            }
            File reAnnotatedCenterVcf = Gatk4MergeSortVcf.sortedVcf.vcf
        }
        File runCsqCenterVcf = select_first([reAnnotatedCenterVcf, runCenterVcf])
        call mergeVcf.PrepParticipantCalls {
            input:
                tumorBarcodeAliquot = tumorBarcodeAliquots[i],
                centerVcf = runCsqCenterVcf,
                referenceFa = referenceFa,
                nonNormalKind = nonNormalKinds[i],
                index = i
        }
    }
    
    # merge into multiNonNormalAliqout and normal VCF
    # annotate with allele read counts
    call mergeAliquots.MergeAliquots {
        input:
            local = local,
            centerVcfs = PrepParticipantCalls.mnvVcf,
            tumorBarcodeAliquots = tumorBarcodeAliquots,
            normalBarcodeAliquot = normalBarcodeAliquot,
            tumorFinalBams = tumorFinalBams,
            normalFinalBam = normalFinalBam,
            participantId = participantId,
            listOfChroms = listOfChroms,
            referenceFa = referenceFa,
            cosmicCensus = cosmicCensus,
            reAnnotate = reAnnotate,
            serviceAccountKey = serviceAccountKey,
            gcpProject = gcpProject
    }
    
    # Make VCFs that are split by pair
    scatter (i in range(length(tumorBarcodeAliquots))) {
    
        call mergeVcf.CompressIndexVcf as intersectionParticipantCompress {
        input:
            vcf = MergeAliquots.intersectionVcf,
            memoryGb = 4
        }
        
    
        Int diskSize = ceil( size(SelectVcfSamples.nonNormalVcf, "GB") * 2) + ceil(size(pairedIntersectionSelectVcfSamples.nonNormalVcf)) + 5
    
    
        String unionPairedVcfPath = "~{tumorBarcodeAliquots[i]}--~{normalBarcodeAliquot}.union.vcf"
        call mergeVcf.SelectVcfSamples as pairedUnionSelectVcfSamples {
            input:
                retainedSamples = [tumorBarcodeAliquots[i], normalBarcodeAliquot],
                centerVcf = MergeAliquots.unionVcf,
                referenceFa = referenceFa,
                nonNormalVcfPath = "~{unionPairedVcfPath}"
        }
        String intersectionPairedVcfPath = "~{tumorBarcodeAliquots[i]}--~{normalBarcodeAliquot}.intersection.vcf"
        call mergeVcf.SelectVcfSamples as pairedIntersectionSelectVcfSamples {
            input:
                retainedSamples = [tumorBarcodeAliquots[i], normalBarcodeAliquot],
                centerVcf = MergeAliquots.intersectionVcf,
                referenceFa = referenceFa,
                nonNormalVcfPath = "~{intersectionPairedVcfPath}"
        }
    }
    
    output {
        Array[File] unionPairedVcf = pairedUnionSelectVcfSamples.nonNormalVcf
        Array[File] intersectionPairedVcf = pairedIntersectionSelectVcfSamples.nonNormalVcf
        File unionVcf = MergeAliquots.unionVcf
        File intersectionVcf = MergeAliquots.intersectionVcf
        File mutect2Genotyped = MergeAliquots.mutect2Genotyped
        File mutect2UnfilteredGenotyped = MergeAliquots.mutect2UnfilteredGenotyped

    }
}