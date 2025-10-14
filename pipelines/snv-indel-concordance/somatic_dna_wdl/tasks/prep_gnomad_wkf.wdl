version 1.0

import "../wdl_structs.wdl"
import "../baf/baf.wdl" as baf
import "utils.wdl" as utils

workflow PrepGnomad {
    input {
        String chrom
        String vcfUri = "gs://gcp-public-data--gnomad/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.sites.chr~{chrom}.vcf.bgz"
        String indexUri = "gs://gcp-public-data--gnomad/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.sites.chr~{chrom}.vcf.bgz.tbi"
        
        String vcfExomeUri = "gs://gcp-public-data--gnomad/release/2.1.1/liftover_grch38/vcf/exomes/gnomad.exomes.r2.1.1.sites.~{chrom}.liftover_grch38.vcf.bgz"
        String vcfExomeIndex = "gs://gcp-public-data--gnomad/release/2.1.1/liftover_grch38/vcf/exomes/gnomad.exomes.r2.1.1.sites.~{chrom}.liftover_grch38.vcf.bgz.tbi"
        Int diskSize = 300
    }
    call utils.DownloadFile as vcfDownloadFile {
        input:
            uri = vcfUri,
            diskSize = diskSize
    }
    call utils.DownloadFile as indexDownloadFile {
        input:
            uri = indexUri,
            diskSize = diskSize
    }
    IndexedVcf referenceVcf = object {
        vcf : vcfDownloadFile.file,
        index : indexDownloadFile.file
    }
    
    # Filter for PASS variants and filter down to AF, AN and nhomalt
    # for non_cancer, XX, XY and full cohorts
    call utils.FilterFile {
        input:
            referenceVcf = referenceVcf,
            filteredReferenceVcfPath = "gnomad.genomes.v3.1.2.pass.sites.chr~{chrom}.vcf.gz",
            diskSize = diskSize
    }
    # Also create biallelic SNP file
    call utils.FilterForBiallelicSnps {
        input:
            gnomadBiallelicPath = "gnomad.genomes.v3.1.2.biallelic.pass.sites.chr~{chrom}.vcf.gz",
            filteredReferenceVcf = FilterFile.filteredReferenceVcf,
            diskSize = 30
    }
    
    # EXOME
    call utils.DownloadFile  as vcfExomeDownloadFile {
        input:
            uri = vcfExomeUri,
            diskSize = diskSize
    }
    call utils.DownloadFile as indexExomeDownloadFile {
        input:
            uri = indexUri,
            diskSize = diskSize
    }
    IndexedVcf referenceExomeVcf = object {
        vcf : vcfExomeDownloadFile.file,
        index : indexExomeDownloadFile.file
    }
    
    # Filter for PASS variants and filter down to AF, AN and nhomalt
    # for non_cancer, XX, XY and full cohorts
    call utils.FilterExomeFile {
        input:
            referenceVcf = referenceExomeVcf,
            filteredReferenceVcfPath = "gnomad.genomes.v3.1.2.pass.sites.chr~{chrom}.vcf.gz",
            diskSize = diskSize
    }
    
    
    output {
        IndexedVcf gnomadGenomes = FilterFile.filteredReferenceVcf
        IndexedVcf gnomadExomes = FilterExomeFile.filteredReferenceVcf
        IndexedVcf gnomadBiallelic = FilterForBiallelicSnps.gnomadBiallelic
    }
}
