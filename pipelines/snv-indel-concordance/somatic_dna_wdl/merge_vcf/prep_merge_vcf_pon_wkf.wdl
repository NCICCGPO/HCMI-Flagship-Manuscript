version 1.0

import "merge_vcf.wdl" as merge_vcf
import "../wdl_structs.wdl"

workflow PrepMergeVcfPon {
    input {
        String tumor
        String tool
        File callerVcf
        IndexedReference referenceFa
    }
    
    call merge_vcf.RenameMetadata {
        input:
            callerVcf = callerVcf,
            tool=tool,
            pairName=tumor
    }
    if (tool == 'manta') {
        call merge_vcf.MergePrepSupport {
            input:
                pairName = tumor,
                tool = tool,
                renameMetaVcf = RenameMetadata.renameMetaVcf
        }
    }
    if (tool != 'manta') {
        call merge_vcf.MergePrep {
            input:
                pairName = tumor,
                tool = tool,
                renameMetaVcf = RenameMetadata.renameMetaVcf
        }

    }
    
    File uncompressedVcf = select_first([MergePrepSupport.prepCallerVcf, MergePrep.prepCallerVcf])
    
    call merge_vcf.RenameVcfPon {
        input:
            tumor = tumor,
            tool = tool,
            prepCallerVcf = uncompressedVcf
    }
    
    call merge_vcf.CompressVcf as renameCompressVcf {
        input:
            vcf = RenameVcfPon.renameVcf,
            memoryGb = 4
    }
    
    call merge_vcf.IndexVcf as renameIndexVcf {
        input:
            vcfCompressed = renameCompressVcf.vcfCompressed
    }
    
    call merge_vcf.SplitMultiAllelic as prepSplitMultiAllelic {
        input:
            pairName = tumor,
            referenceFa = referenceFa,
            vcfCompressedIndexed = renameIndexVcf.vcfCompressedIndexed,
            splitVcfPath = sub(basename(renameIndexVcf.vcfCompressedIndexed.vcf, ".gz"), ".vcf$", ".split_ma.vcf")
    }
    
    call merge_vcf.SplitMnv {
        input:
            tool = tool,
            mnvVcfPath = sub(basename(renameIndexVcf.vcfCompressedIndexed.vcf, ".gz"), ".vcf$", ".split_mnv.vcf"),
            splitVcf = prepSplitMultiAllelic.splitVcf
            
    }
    
    if (tool == 'svaba') {
        call merge_vcf.RemoveContig {
            input:
                mnvVcfPath = sub(basename(renameIndexVcf.vcfCompressedIndexed.vcf, ".gz"), ".vcf$", ".split_no_contig.vcf"),
                removeChromVcf = SplitMnv.mnvVcf
        }
    }
    
    call merge_vcf.Gatk4MergeSortVcf {
        input:
            tempVcfs = [select_first([RemoveContig.removeContigVcf, SplitMnv.mnvVcf])],
            sortedVcfPath = sub(basename(select_first([RemoveContig.removeContigVcf, SplitMnv.mnvVcf])), "$", ".gz"),
            referenceFa = referenceFa,
            gzipped = true,
            threads = 4,
            memoryGb = 8,
            diskSize = 20
            
    }
    
    output {
        IndexedVcf preppedVcf = Gatk4MergeSortVcf.sortedVcf
    }
}
