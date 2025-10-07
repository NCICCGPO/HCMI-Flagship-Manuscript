#!/bin/bash

# script to take genotype file as called from pileup and format to be used by PLINK
# then to be merged w 1000Genomes references
# for eventual use with the ADMIXTURE software

# takes as input:
# full file path and name to the sample called genotypes
#    assumes file includes 6 cols: chr, pos, rsID, alleleA, alleleB, genotype formatted as 0/0,0/1,1/0,1/1 or NA for missing
# output file path/name

uuid=$(uuidgen)

INPUTFILE=$1
OUTFILE=$2
MARKERSET=$OUTFILE"_markerSet"
HEADERFILE=$3
HEADERFILE2=$4

# first, get exact marker list
awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5}' <$INPUTFILE >$MARKERSET

cat $HEADERFILE $MARKERSET > $OUTFILE

cut -f6 $INPUTFILE > $OUTFILE.tmpFile.txt
tmpName=$(basename $INPUTFILE)
SNAME=${tmpName/.genotypes.txt/}
echo $SNAME | cat - $OUTFILE.tmpFile.txt > $OUTFILE.tmpFile2.txt
paste -d "\t" $OUTFILE $OUTFILE.tmpFile2.txt > $OUTFILE.tmpFile.txt
mv $OUTFILE.tmpFile.txt $OUTFILE

sed 's/0\/0/2/g' $OUTFILE > $OUTFILE.tmpFile.txt
sed 's/0\/1/1/g' $OUTFILE.tmpFile.txt > $OUTFILE
sed 's/1\/0/1/g' $OUTFILE > $OUTFILE.tmpFile.txt
sed 's/1\/1/0/g' $OUTFILE.tmpFile.txt > $OUTFILE

awk '{print $1"\t"$3"\t"0"\t"$2"\t"$4"\t"$5"\t"$6}' <$OUTFILE >$OUTFILE.tmpFile.txt
mv $OUTFILE.tmpFile.txt $OUTFILE
awk 'NR>1' $OUTFILE >$OUTFILE"_noHdr"
echo "0_"$SNAME | paste -d "\t" $HEADERFILE2 - > ${HEADERFILE2}_${uuid}".tmp"
cat ${HEADERFILE2}_${uuid}.tmp $OUTFILE"_noHdr" > $OUTFILE

rm $OUTFILE"_noHdr"
rm ${HEADERFILE2}_${uuid}.tmp
rm $OUTFILE.tmpFile2.txt
rm $MARKERSET
