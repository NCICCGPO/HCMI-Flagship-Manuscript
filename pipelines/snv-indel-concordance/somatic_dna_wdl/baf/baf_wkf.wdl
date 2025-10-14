version 1.0

import "../wdl_structs.wdl"
import "baf.wdl"
import "../test/tests.wdl" as tests
import "../tasks/utils.wdl" as utils
import "../merge_vcf/merge_vcf.wdl" as mergeVcf
import "../calling/calling.wdl" as calling


workflow Baf {
    # command 
    input {
        Boolean local = false
        String sampleId
        String pairName
        Bam normalFinalBam
        Bam tumorFinalBam
        File? germlineVcf
        IndexedReference referenceFa
        
        File? serviceAccountKey
        String? gcpProject
    }
    if (!defined(germlineVcf)) {
        call utils.CreateBlankFile {
            input:
                fileId = "~{sampleId}.fake_filler_"
        }
    }
    File inputGermlineVcf = select_first([germlineVcf, CreateBlankFile.blankFile])
    # gather split VCF filenames to avoid glob that can fail on prem
    Array[String] suffixes = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29", "30", "31", "32", "33", "34", "35", "36", "37", "38", "39", "40", "41", "42", "43", "44", "45", "46", "47", "48", "49", "50"]
    if (!local) {
        Int cloudMaxSplits = 30
        Int cloudMinSplits = 2
    }
    if (local) {
        Int hpcMaxSplits = 2
        Int hpcMinSplits = 1
    }
    String prefix = "~{sampleId}.haplotypecaller.gatk"
    String additionalSuffix = ".vcf"
    Array[String] mnvVcfPaths = ["~{prefix}.mnv~{additionalSuffix}"]
    scatter (index in range(length(suffixes))) {
        String suffix = suffixes[index]
        String splitChromVcfPaths = "~{prefix}.~{suffix}~{additionalSuffix}"
    }
    Array[String] splitVcfPaths = flatten([mnvVcfPaths, splitChromVcfPaths])
    if (defined(germlineVcf)) {
        call utils.SplitVcf {
            input:
                vcf=inputGermlineVcf,
                prefix=prefix,
                diskSize = (ceil(size(inputGermlineVcf, "GB")) * 3) + 10,
                maxRows = 1000,
                minSplits = select_first([cloudMinSplits, hpcMinSplits]),
                maxSplits = select_first([cloudMaxSplits, hpcMaxSplits]),
                splitVcfPaths = splitVcfPaths
        }
        Array[File] splitVcfs = select_all(SplitVcf.splitVcfs)
        scatter(vcf in splitVcfs) {
            call baf.FilterForHetSnps {
                input:
                    sampleId = sampleId,
                    referenceFa = referenceFa,
                    germlineVcf = vcf
            }
            
            call baf.FilterBaf {
                input:
                    sampleId = sampleId,
                    hetVcf = FilterForHetSnps.hetVcf
            }
            
            call mergeVcf.VcfToBed {
                input:
                    pairName = pairName,
                    chrom = 'split',
                    candidateChromVcf = FilterBaf.knownHetVcf,
                    candidateChromBedPath = "~{pairName}.FilterBaf.bed"
            }
            
            call calling.SliceBam as sliceBamTumor {
                input:
                    chrom = 'split',
                    sampleId = pairName,
                    chromBed = VcfToBed.candidateChromBed,
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
                    sampleId = sampleId,
                    chromBed = VcfToBed.candidateChromBed,
                    finalBam = normalFinalBam,
                    threads = 4,
                    memoryGb = 4,
                    diskSize = 10,
                    serviceAccountKey = serviceAccountKey,
                    gcpProject = gcpProject
            }
    
            Int diskSize = ceil( size(sliceBamNormal.lancetBamSlice.bam, "GB") + size(sliceBamTumor.lancetBamSlice.bam, "GB")) + 10
            call baf.AlleleCountsLocalize {
                input:
                    pairName = pairName,
                    referenceFa = referenceFa,
                    normalFinalBam = sliceBamNormal.lancetBamSlice,
                    tumorFinalBam = sliceBamTumor.lancetBamSlice,
                    knownHetVcf = FilterBaf.knownHetVcf,
                    diskSize = diskSize
            }
            call baf.CalcBaf {
                input:
                    pairName = pairName,
                    alleleCountsTxt = AlleleCountsLocalize.alleleCountsTxt
            }
        }
        
        call tests.ConcateTables {
            input:
                tables = CalcBaf.bafTxt,
                outputTablePath = "~{pairName}.haplotypecaller.gatk.baf.txt"
        }
    }
    
    output {
        File? alleleCountsTxt = ConcateTables.outputTable
    }
}
