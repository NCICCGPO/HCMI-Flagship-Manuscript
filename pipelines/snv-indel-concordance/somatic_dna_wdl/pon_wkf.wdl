version 1.0

import "calling/calling.wdl" as callingTasks
import "calling/mutect2_pon_wkf.wdl" as mutect2Pon
import "calling/manta_pon_wkf.wdl" as mantaPon
import "merge_vcf/merge_vcf_pon_wkf.wdl" as MergerVcfPonWgs
import "merge_vcf/merge_vcf_exome_pon_wkf.wdl" as MergerVcfPonExome
import "annotate/annotate.wdl" as annotate
import "annotate/variantEffectPredictor.wdl" as variantEffectPredictor
import "merge_vcf/merge_vcf.wdl" as mergeVcf


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

import "wdl_structs.wdl"

workflow CallingPon {
    # command
    #   Call variants in BAMs
    #   merge and filter raw VCFs
    #   annotate
    input {
        Boolean production = true
        String library
        Array[SampleBamInfo]+ tumorInfos
        # strelka2
        File strelkaJsonLog
        File configureStrelkaSomaticWorkflow
        #   mutect2
        File mutectJsonLog
        Array[String]+ listOfChroms
        Array[String]+ callerIntervals
        File invertedIntervalListBed
        IndexedReference referenceFa
        #   Manta
        IndexedTable callRegions
        File mantaJsonLog
        File mutectJsonLogFilter
        File refCache
        Boolean highMem = false
        
        # annotation:
        IndexedVcf cosmicCoding
        IndexedVcf cosmicNoncoding
        
        # Public
        File cancerResistanceMutations
        File vepCache
        File annotations
        File plugins
        String vepGenomeBuild
        IndexedReference vepFastaReference

        # NYGC-only
        IndexedVcf hgmdGene
        IndexedVcf hgmdUd10
        IndexedVcf hgmdPro
        IndexedVcf omimVcf

        # Public
        IndexedVcf chdGenesVcf
        IndexedVcf chdEvolvingGenesVcf
        IndexedVcf chdWhitelistVcf
        IndexedVcf deepIntronicsVcf
        IndexedVcf clinvarIntronicsVcf
        IndexedVcf masterMind

        
        Int vepDiskSize = ceil(size(vepCache, "GB") + size(plugins, "GB") + size(annotations, "GB") + size(chdGenesVcf.vcf, "GB") + size(chdEvolvingGenesVcf.vcf, "GB") + size(chdWhitelistVcf.vcf, "GB") + size(deepIntronicsVcf.vcf, "GB") + size(clinvarIntronicsVcf.vcf, "GB") + size(masterMind.vcf, "GB")) + 200
        
    }
    
    scatter(tumorInfo in tumorInfos) {
        
        call mutect2Pon.Mutect2Pon {
            input:
                mutectJsonLog = mutectJsonLog,
                mutectJsonLogFilter = mutectJsonLogFilter,
                tumor = tumorInfo.sampleId,
                callerIntervals = callerIntervals,
                referenceFa = referenceFa,
                tumorFinalBam = tumorInfo.finalBam,
                highMem = highMem,
                invertedIntervalListBed = invertedIntervalListBed,
                library = library
        }
        
        if (library == 'WGS') {
            call mantaPon.MantaPon {
                input:
                    mantaJsonLog = mantaJsonLog,
                    tumor = tumorInfo.sampleId,
                    callRegions = callRegions,
                    referenceFa = referenceFa,
                    tumorFinalBam = tumorInfo.finalBam,
                    highMem = highMem
            }
            
            call MergerVcfPonWgs.MergeVcfPonWgs {
                        input:
                            tumor = tumorInfo.sampleId,
                            filteredMantaSV = MantaPon.filteredMantaSV,
                            mutect2 = Mutect2Pon.mutect2,
                            
                            referenceFa = referenceFa,
                            listOfChroms = listOfChroms
            }
        }
        
        if (library == 'Exome') {
            call MergerVcfPonExome.MergeVcfPonExome {
                        input:
                            tumor = tumorInfo.sampleId,
                            mutect2 = Mutect2Pon.mutect2,
                            
                            referenceFa = referenceFa,
                            listOfChroms = listOfChroms
            }
        }
        
        File mergedVcf = select_first([MergeVcfPonWgs.mergedVcf, MergeVcfPonExome.mergedVcf])
        
        call mergeVcf.CompressVcf as unannotatedCompressVcf {
            input:
                vcf = mergedVcf
        }
    
        call mergeVcf.IndexVcf as unannotatedIndexVcf {
            input:
                vcfCompressed = unannotatedCompressVcf.vcfCompressed
        }
    
        if (production) {
            call variantEffectPredictor.vepPublicSvnIndel as productionVepSvnIndel {
                input:
                    pairName = tumorInfo.sampleId,
                    unannotatedVcf = unannotatedIndexVcf.vcfCompressedIndexed,
                    vepCache = vepCache,
                    annotations = annotations,
                    plugins = plugins,
                    vepGenomeBuild = vepGenomeBuild,
                    vepFastaReference = vepFastaReference,
                    # Public
                    chdGenesVcf = chdGenesVcf,
                    chdEvolvingGenesVcf = chdEvolvingGenesVcf,
                    chdWhitelistVcf = chdWhitelistVcf,
                    deepIntronicsVcf = deepIntronicsVcf,
                    clinvarIntronicsVcf = clinvarIntronicsVcf,
                    masterMind = masterMind,
                    cosmicCoding = cosmicCoding,
                    cosmicNoncoding = cosmicNoncoding,
                    diskSize = vepDiskSize,
                    memoryGb = 8
            }
        }
        
        if (!production) {
            if ( size(omimVcf.vcf) > 0 ) {
                call variantEffectPredictor.vepSvnIndel as notProductionVepSvnIndel{
                    input:
                        pairName = tumorInfo.sampleId,
                        unannotatedVcf = unannotatedIndexVcf.vcfCompressedIndexed,
                        vepCache = vepCache,
                        annotations = annotations,
                        plugins = plugins,
                        vepGenomeBuild = vepGenomeBuild,
                        vepFastaReference = vepFastaReference,
                        # NYGC-only
                        hgmdGene = hgmdGene,
                        hgmdUd10 = hgmdUd10,
                        hgmdPro = hgmdPro,
                        omimVcf = omimVcf,
                        # Public
                        chdGenesVcf = chdGenesVcf,
                        chdEvolvingGenesVcf = chdEvolvingGenesVcf,
                        chdWhitelistVcf = chdWhitelistVcf,
                        deepIntronicsVcf = deepIntronicsVcf,
                        clinvarIntronicsVcf = clinvarIntronicsVcf,
                        masterMind = masterMind,
                        cosmicCoding = cosmicCoding,
                        cosmicNoncoding = cosmicNoncoding,
                        diskSize = vepDiskSize,
                        memoryGb = 8
                }
            }
        }
    
        File finalVcfAnnotatedVep = select_first([notProductionVepSvnIndel.vcfAnnotatedVep, productionVepSvnIndel.vcfAnnotatedVep])
        
    }
 
    output {
        # temp
        # General
        Array[File] vcfAnnotatedVep = finalVcfAnnotatedVep
        
        # Mutect2
        Array[File] mutect2 = Mutect2Pon.mutect2
        # Manta 
        Array[File?] filteredMantaSV = MantaPon.filteredMantaSV

    }
}
