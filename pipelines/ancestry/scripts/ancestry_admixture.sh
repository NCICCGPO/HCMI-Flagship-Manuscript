#!/bin/bash
#SBATCH --mem=80G
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=1

# set -euf -o pipefail

script_root=`dirname $(readlink -f "$0")`
USAGE_MESSAGE="Usage: [ -bam <BAM> | -gvcf <GVCF> | -gt <GENOTYPE> ] -sample <sample_id> -outdir <output_dir> -seq_type <wgs|exome|impact|rna> [ -subpops ]"
subpops=0

while true; do
    case $1 in
    -h|-\?|--help)
        echo -e $USAGE_MESSAGE >&2
        exit 0
        ;;
    -bam)
        shift
        bam=$1
        ;;
    -gvcf)
        shift
        gvcf=$1
        ;;
    -gt)
        shift
        gt=$1
        ;;
    -sample)
        shift
        sample=$1
        ;;
    -outdir)
        shift
        outdir=$1
        ;;
    -subpops)
        shift
        subpops=1
        ;;
    -?*)
        printf 'WARNING: Unknown option (ignored): %s\n' "$1" >&2
        echo -e $USAGE_MESSAGE >&2
        ;;
    *)
        break
    esac
    shift
done


if [[ -z "$outdir" ]]; then
    echo Please specify output directory
    echo -e $USAGE_MESSAGE >&2
    exit -1
fi

ref_ver=b38
seq_type=wgs
tmpdir=${outdir}/tmp_ancestry_cm_${sample}
mkdir -p ${tmpdir}

if [[ "${ref_ver}" == "b37" ]] || [[ "${ref_ver}" == "b38" ]]; then
    echo "Supported reference version."
    dbsnp=`jq .${ref_ver}.dbSNP ${script_root}/../resources/cm/markers/path_to_reference_and_markers_for_genotyping.json`
    markerprefix=`jq .${ref_ver}.${seq_type}.markerFilePrefix ${script_root}/../resources/cm/markers/path_to_reference_and_markers_for_genotyping.json`
    dbsnp=`echo $dbsnp | sed 's/"//g'`
    markerprefix=`echo $markerprefix | sed 's/"//g'`
    echo $dbsnp
    markers=${markerprefix}.prune.in
    markerstxt=${markerprefix}.markers.txt
    kgp_plink=$markerprefix
    if [[ "${ref_ver}" == "b37" ]] && [[ "${seq_type}" == "rna" ]]; then ## Comment from Kanika on 11/26/2019. Very hacky way of running on RNA at the moment. Needs to be fixed.
        seq_type="exome"
        dbsnp=`jq .${ref_ver}.dbSNP ${script_root}/../resources/cm/markers/path_to_reference_and_markers_for_genotyping.json`
        markerprefix=`jq .${ref_ver}.${seq_type}.markerFilePrefix ${script_root}/../resources/cm/markers/path_to_reference_and_markers_for_genotyping.json`
        dbsnp=`echo $dbsnp | sed 's/"//g'`
        markerprefix=`echo $markerprefix | sed 's/"//g'`
        markers=${markerprefix}.prune.in
        markerstxt=${markerprefix}.markers.txt
        kgp_plink=$markerprefix
        seq_type="rna"
    fi  
else
    echo "Unknown reference version."
    exit 1
fi

if [ -z "$bam" ] && [ -z "$gvcf" ] && [ -z "$gt" ]; then
    echo "Either bam file or gVCF file required."
    echo "Usage: [ -bam <BAM> | -gvcf <GVCF> | -gt <GENOTYPE_FILE> ] -sample <sample_id> -ref_ver <b38|b37> -outdir <output_dir> -seq_type <wgs|exome>"
    exit -1
fi

if [ "${seq_type}" == "wgs" ] || [ "${seq_type}" == "exome" ] || [ "${seq_type}" == "impact" ] || [ "${seq_type}" == "rna" ]; then
    echo "Supported sequencing type."
else
    echo "Unsupported sequencing type. Please choose either 'wgs' or 'exome' or 'impact' or 'rna'."
    exit 1
fi

if [[ ! -z "$bam" ]]; then
    python ${script_root}/../resources/cm/scripts/run_gatk_pileup_for_sample.py \
        -B ${bam} \
        -O ${tmpdir}/${sample}.pileup \
        -G ${ref_ver} \
        -T ${seq_type}

    python ${script_root}/../resources/cm/scripts/calculate_genotypes_for_sample.py \
        -P ${tmpdir}/${sample}.pileup \
        -O ${tmpdir}/${sample}.genotypes.txt \
        -G ${ref_ver} \
        -T ${seq_type}


    ${script_root}/../resources/cm/scripts/batch_parseGenos.sh \
        ${tmpdir}/${sample}.genotypes.txt \
        ${tmpdir}/${sample}.toGenos \
        ${script_root}/../resources/header.txt \
        ${script_root}/../resources/header2.txt
elif [[ ! -z "$gvcf" ]]; then
    python ${script_root}/../resources/cm/scripts/gvcf_to_gt.py \
        -i ${gvcf} \
        -m ${markerstxt} \
        -d ${dbsnp} \
        -r ${ref_ver} \
        -o ${tmpdir}/${sample}.toGenos
else
   cp $gt ${tmpdir}/${sample}.genotypes.txt
   ${script_root}/../resources/cm/scripts/batch_parseGenos.sh \
        ${tmpdir}/${sample}.genotypes.txt \
        ${tmpdir}/${sample}.toGenos
fi

if [[ "${ref_ver}" == "b37" ]] && [[ "${seq_type}" == "rna" ]]; then ## Needs to be improved
    mv ${tmpdir}/${sample}.toGenos ${tmpdir}/${sample}.with_chr_prefix.toGenos
    sed 's/^chr//' ${tmpdir}/${sample}.with_chr_prefix.toGenos > ${tmpdir}/${sample}.toGenos
fi

echo BAM/VCF to GT done
kgp_plink=`echo ${kgp_plink} | sed 's/"//g'`
${script_root}/../resources/cm/scripts/batch_mergeWRef.sh \
    -genos ${tmpdir}/${sample}.toGenos \
    -outfile ${tmpdir}/${sample}.plink \
    -kgp_plink ${kgp_plink} 
echo batch_mergeWRef.sh done

${script_root}/../resources/cm/scripts/batch_runADMIXTURE.sh \
    -script_root ${script_root} \
    -plink ${tmpdir}/${sample}.plink_w1kgRef
echo batch_runADMIXTURE.sh done

numpops=`grep -vw "-" ${tmpdir}/${sample}.plink_w1kgRef.pop | sort | uniq | wc -l`
cp ${tmpdir}/${sample}.plink_w1kgRef.admixtureRes.${numpops}pops.studyOnly ${outdir}/.

## Plot PCA ##
RSCRIPT=${script_root}/../resources/cm/scripts/plotPCA.R

plink --bfile ${tmpdir}/${sample}.plink_w1kgRef \
    --pca \
    --out ${tmpdir}/${sample}.plink_w1kgRef

R --vanilla --args ${tmpdir}/${sample}.plink_w1kgRef.eigenvec ${tmpdir}/${sample}.plink_w1kgRef.admixtureRes.${numpops}pops ${tmpdir}/${sample}.plink_w1kgRef <$RSCRIPT
cp ${tmpdir}/${sample}.plink_w1kgRef.pca.pdf ${outdir}/.

if [ $subpops == 1 ]; then
${script_root}/../resources/cm/scripts/batch_runADMIXTURE.subpops.sh \
    -script_root ${script_root} \
    -plink ${tmpdir}/${sample}.plink_w1kgRef
cp ${tmpdir}/${sample}.plink_w1kgRef.admixtureRes.subpops.studyOnly ${outdir}/.
fi
#rm -r ${tmpdir}
