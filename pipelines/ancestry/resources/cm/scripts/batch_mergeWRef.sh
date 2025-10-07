#!/bin/bash

# script to merge genotype file of study samples with 1kg reference file
# takes as input:
# genotype file of study samples
# assumes ONLY ONE study sample in genotype file

# assumes genotypes are in format 0/0,0/1,1/0,1/1 or NA for missing
# assumes merging w 1kg PLINK file


while [[ $# -gt 1 ]]; do
    key=${1}
    case ${key} in
        -script_root)
            script_root=${2};;
        -genos)
            genos=${2};;
        -outfile)
            outfile=${2};;
        -kgp_plink)
            kgp_plink=${2};;
        *)
            echo "Unknown option: ${key}"
            echo "Usage: -genos <sample>.toGenos -ref_ver <b37|b38> -outfile <outfile> "
            exit 65
        ;;
    esac
    shift
    shift
done

# Don't ask me why there is plink1 and 2. It always existed.
REF=${kgp_plink}

# read study genotypes into plink2 'traw' format
# need to make a .fam file first

SID=$(awk 'NR==1 {print $7}' $genos)
NSID=$(echo "${SID:2:${#SID}-2}")
awk -v var="$NSID" 'NR==1 {print 0"\t"var"\t"0"\t"0"\t"0"\t-9"}' <$genos >$outfile.inputfam

plink --import-dosage $genos skip0=1 skip1=2 skip2=0 chr-col-num=1 pos-col-num=4 format=1 id-delim="_" \
	--fam $outfile.inputfam \
	--no-pheno \
	--make-bed \
	--geno 0.98 \
	--set-missing-var-ids @:# \
	--out $outfile

plink --bfile $outfile \
       --bmerge $REF \
       --out $outfile"_w1kgRef"

if [ -f $outfile"_w1kgRef.missnp" ]; then
    plink --bfile $outfile \
        --exclude $outfile"_w1kgRef.missnp" \
        --make-bed \
        --out $outfile"_wExcl"
else
    plink --bfile $outfile \
        --make-bed \
        --out $outfile"_wExcl"
fi

# subset the 1kgRef to the variants included in the $outfile_wExcl
awk '{print $2}' <$outfile"_wExcl.bim" >$outfile"_vars"

plink --bfile $REF \
       --extract $outfile"_vars" \
       --make-bed \
       --out $outfile"_refVars"

plink --bfile $outfile"_wExcl" \
       --bmerge $outfile"_refVars" \
       --out $outfile"_w1kgRef"

plink --bfile $outfile"_w1kgRef" \
       --geno 0.98 \
       --make-bed \
       --out $outfile"_w1kgRef"
