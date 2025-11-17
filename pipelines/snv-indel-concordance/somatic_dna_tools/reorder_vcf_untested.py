from pysam import VariantFile
import sys


def reorder_vcf(vcf_in, vcf_out, normal, tumor):
    """
    Reorders VCF columns to place 'normal' and 'tumor' samples in the desired order.

    Args:
        vcf_in: Path to the input VCF file.
        vcf_out: Path to the output VCF file.
        normal: Sample name for the normal sample.
        tumor: Sample name for the tumor sample.
    """

    with VariantFile(vcf_in) as vcf_reader, VariantFile(vcf_out, 'w', header=vcf_reader.header) as vcf_writer:
        # Get sample indices for normal and tumor
        normal_idx = vcf_reader.header.index(normal)
        tumor_idx = vcf_reader.header.index(tumor)

        # Check if normal and tumor samples exist
        if normal_idx == -1 or tumor_idx == -1:
            raise ValueError(f"Samples '{normal}' and '{tumor}' not found in VCF header.")

        # Reorder header and samples if necessary
        new_header = vcf_reader.header.copy()
        new_header.samples = [normal, tumor]
        vcf_writer.header = new_header

        # Write records with reordered samples
        for record in vcf_reader:
            new_genotypes = [record.samples[normal_idx], record.samples[tumor_idx]]
            record.samples = new_genotypes
            vcf_writer.write(record)


if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: reorder_vcf.py VCF_IN VCF_OUT NORMAL TUMOR")
        sys.exit(1)

    vcf_in, vcf_out, normal, tumor = sys.argv[1:]
    reorder_vcf(vcf_in, vcf_out, normal, tumor)
