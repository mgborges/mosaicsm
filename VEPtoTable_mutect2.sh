function VEPtoTable
{
VCF=$1
VEP=$VCF.VEP
grep -v ^# $VEP > $VCF.CABECALHO

#java -jar /home/murilo/snpEff/SnpSift.jar extractFields $VCF CHROM POS ALT REF "GEN[*].GT" "GEN[*].AD" -e "." -s ";" > $VCF.tmp
java -jar /home/murilo/snpEff/SnpSift.jar extractFields $VCF CHROM POS ALT REF "GEN[*].GT" "GEN[*].AD" -e "." -s ";" > $VCF.tmp

while read -r line
do
IMPACT=`echo $line | awk -F "IMPACT=" '{print $2}' | awk -F ";" '{print $1}'`
if [ "$IMPACT" == "" ]
then
IMPACT="."
fi
SYMBOL=`echo $line | awk -F "SYMBOL=" '{print $2}' | awk -F ";" '{print $1}'`
if [ "$SYMBOL" == "" ]
then
SYMBOL="."
fi
BIOTYPE=`echo $line | awk -F "BIOTYPE=" '{print $2}' | awk -F ";" '{print $1}'`
if [ "$BIOTYPE" == "" ]
then
BIOTYPE="."
fi
GMAF=`echo $line | awk -F "GMAF=" '{print $2}' | awk -F ";" '{print $1}'`
if [ "$GMAF" == "" ]
then
GMAF="."
fi
SIFT=`echo $line | awk -F "SIFT=" '{print $2}' | awk -F ";" '{print $1}'`
if [ "$SIFT" == "" ]
then
SIFT="."
fi
PolyPhen=`echo $line | awk -F "PolyPhen=" '{print $2}' | awk -F ";" '{print $1}'`
if [ "$PolyPhen" == "" ]
then
PolyPhen="."
fi
otherinfo=`echo $line | cut -f 1-13 -d " "`

CHR=`echo $line | cut -f 2 -d " " | cut -f 1 -d ":"`
POS=`echo $line | cut -f 2 -d " " | cut -f 2 -d ":" | cut -f 1 -d "-"`
ALL=`echo $line | cut -f 3 -d " "`
#CHR="chr$CHR"

verificator=`echo $line | cut -f 2 -d " " | cut -f 2 -d ":" | grep -c "-"`

if [ "$verificator" != "0" ] && [ "$ALL" != "-" ]
then
GT=`grep -w -e ^$CHR$'\t'$POS $VCF.tmp | cut -f 5`
DP=`grep -w -e ^$CHR$'\t'$POS $VCF.tmp | cut -f 6`
elif [ "$verificator" != "0" ] && [ "$ALL" == "-" ]
then
POS=`bc <<< "$POS - 1"`
GT=`grep -w -e ^$CHR$'\t'$POS $VCF.tmp | cut -f 5`
DP=`grep -w -e ^$CHR$'\t'$POS $VCF.tmp | cut -f 6`
elif [ "$verificator" == "0" ] && [ "$ALL" == "-" ]
then
POS=`bc <<< "$POS - 1"`
GT=`grep -w -e ^$CHR$'\t'$POS $VCF.tmp | cut -f 5`
DP=`grep -w -e ^$CHR$'\t'$POS $VCF.tmp | cut -f 6`
else
GT=`grep -w -e ^$CHR$'\t'$POS $VCF.tmp | cut -f 5`
DP=`grep -w -e ^$CHR$'\t'$POS $VCF.tmp | cut -f 6`
fi

REF_TECIDO=`echo $DP | cut -f 2 -d ';' | cut -f 1 -d ','`

calc=`echo $DP | cut -f 2 -d ';' | cut -f 2-10 -d , | sed 's/,/+/g'`


ALT_SANGUE=0

MOSAICO_TECIDO=`grep -w -e ^$CHR$'\t'$POS $VCF | cut -f 11 | cut -f 3 -d :`
MOSAICO_SANGUE=0

BAM_TECIDO=`echo $vcf | cut -f 2 -d - | cut -f 1 -d .`
BAM_SANGUE=`echo $vcf | cut -f 1 -d -`

cobertura_TECIDO=`samtools depth -r $CHR:$POS-$POS ../bams/$BAM_TECIDO*.bam | cut -f 3`
cobertura_SANGUE=`samtools depth -r $CHR:$POS-$POS ../bams/$BAM_SANGUE*.bam | cut -f 3`

PON=`grep -w -e ^$CHR$'\t'$POS $VCF| grep IN_PON -c`

echo $otherinfo $IMPACT $SYMBOL $BIOTYPE $GMAF $SIFT $PolyPhen $PON $GT $DP $MOSAICO_SANGUE $cobertura_SANGUE $MOSAICO_TECIDO $cobertura_TECIDO 
done < $VCF.CABECALHO
rm $VCF.tmp $VCF.CABECALHO
}

LC_ALL=C locale

for vcf in `ls *.vcf`
do
	# VEP $vcf
	echo PROCESSING $vcf ... &&
	VEPtoTable $vcf > $vcf.txt 2> $vcf.txt.log &&
	cat ~murilo/Dropbox/sh_scripts/VEPtoTable_mutect2.header $vcf.txt > $vcf.txt2 &&
	mv  $vcf.txt2 $vcf.txt &&
	echo CONCLUDED WITH $vcf
done

