version 1.0

import "calling.wdl" as calling
import "../wdl_structs.wdl"

workflow Strelka2 {
    # command
    #   run Strelka2 caller
    input {
        String library
        String tumor
        String normal
        String intHVmem = "unlimited"
        IndexedReference referenceFa
        # candidateSmallIndels may need to be a dummy file (if you want to require it for the task)
        IndexedVcf candidateSmallIndels
        File intervalListBed
        Bam normalFinalBam
        IndexedTable callRegions
        Bam tumorFinalBam
        String pairName
        String strelka2SnvPath = "~{pairName}.snv.strelka2.vcf"
        String strelka2IndelPath = "~{pairName}.indel.strelka2.vcf"
        # resources
        Int diskSize = ceil( size(tumorFinalBam.bam, "GB") + size(normalFinalBam.bam, "GB")) + 20
        Int memoryGb = 4
        Int threads = 8
        File configureStrelkaSomaticWorkflow
        File strelkaJsonLog
    }

    if (library == 'WGS') {
        call calling.Strelka2 {
            input:
                configureStrelkaSomaticWorkflow = configureStrelkaSomaticWorkflow,
                intHVmem = intHVmem,
                referenceFa = referenceFa,
                callRegions = callRegions,
                pairName = pairName,
                candidateSmallIndels = candidateSmallIndels,
                normalFinalBam = normalFinalBam,
                tumorFinalBam = tumorFinalBam,
                memoryGb = memoryGb,
                diskSize = diskSize,
                threads = threads
        }
    }
    
    if (library == 'Exome') {
        call calling.Strelka2Exome {
            input:
                configureStrelkaSomaticWorkflow = configureStrelkaSomaticWorkflow,
                intHVmem = intHVmem,
                referenceFa = referenceFa,
                callRegions = callRegions,
                pairName = pairName,
                normalFinalBam = normalFinalBam,
                tumorFinalBam = tumorFinalBam,
                memoryGb = memoryGb,
                diskSize = diskSize,
                threads = threads
        }
        
        call calling.SelectVariants as selectVariantsIndel {
            input:
                referenceFa = referenceFa,
                pairName = pairName,
                vcf =  Strelka2Exome.strelka2Indels,
                intervalListBed = intervalListBed,
                memoryGb = memoryGb,
                diskSize = diskSize,
                threads = threads
        }
        
        call calling.SelectVariants as selectVariantsSnv {
            input:
                referenceFa = referenceFa,
                pairName = pairName,
                vcf =  Strelka2Exome.strelka2Snvs,
                intervalListBed = intervalListBed,
                memoryGb = memoryGb,
                diskSize = diskSize,
                threads = threads
        }
    }

    IndexedVcf strelka2SnvsFinal = select_first([Strelka2.strelka2Snvs, Strelka2Exome.strelka2Snvs])
    IndexedVcf strelka2IndelsFinal = select_first([Strelka2.strelka2Indels, Strelka2Exome.strelka2Indels])
    File strelka2SnvsVcf = select_first([strelka2SnvsFinal.vcf, selectVariantsSnv.outVcf])
    File strelka2IndelsVcf = select_first([strelka2IndelsFinal.vcf, selectVariantsIndel.outVcf])
    
    call calling.AddVcfCommand as strelka2SnvAddVcfCommand {
        input:
            inVcf = strelka2SnvsVcf,
            jsonLog = strelkaJsonLog,
            memoryGb = 4,
            diskSize = 10
    }

    call calling.ReorderVcfColumns as strelka2SnvReorderVcfColumns {
        input:
            tumor = tumor,
            normal = normal,
            rawVcf = strelka2SnvAddVcfCommand.outVcf,
            orderedVcfPath = strelka2SnvPath,
            memoryGb = 4,
            diskSize = 10
    }

    call calling.AddVcfCommand as strelka2IndelAddVcfCommand {
        input:
            inVcf = strelka2IndelsVcf,
            jsonLog = strelkaJsonLog,
            memoryGb = 4,
            diskSize = 10
    }

    call calling.ReorderVcfColumns as strelka2IndelReorderVcfColumns {
        input:
            tumor = tumor,
            normal = normal,
            rawVcf = strelka2IndelAddVcfCommand.outVcf,
            orderedVcfPath = strelka2IndelPath,
            memoryGb = 4,
            diskSize = 10
    }

    output {
        IndexedVcf strelka2Snvs = strelka2SnvsFinal
        IndexedVcf strelka2Indels = strelka2IndelsFinal
        File strelka2Snv = strelka2SnvReorderVcfColumns.orderedVcf
        File strelka2Indel = strelka2IndelReorderVcfColumns.orderedVcf
    }
}
