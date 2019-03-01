function VEPtoTable
{
VCF=$1
VEP=$VCF.VEP
grep -v ^# $VEP > $VCF.CABECALHO

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

######################################################################################

if [ "$verificator" != "0" ] && [ "$ALL" != "-" ]
then
GT_SANGUE=`grep -w -e ^$CHR$'\t'$POS $VCF | cut -f 10 | cut -f 1 -d :`
GT_TECIDO=`grep -w -e ^$CHR$'\t'$POS $VCF | cut -f 11 | cut -f 1 -d :`
elif [ "$verificator" != "0" ] && [ "$ALL" == "-" ]
then
POS=`bc <<< "$POS - 1"`
GT_SANGUE=`grep -w -e ^$CHR$'\t'$POS $VCF | cut -f 10 | cut -f 1 -d :`
GT_TECIDO=`grep -w -e ^$CHR$'\t'$POS $VCF | cut -f 11 | cut -f 1 -d :`
elif [ "$verificator" == "0" ] && [ "$ALL" == "-" ]
then
POS=`bc <<< "$POS - 1"`
GT_SANGUE=`grep -w -e ^$CHR$'\t'$POS $VCF | cut -f 10 | cut -f 1 -d :`
GT_TECIDO=`grep -w -e ^$CHR$'\t'$POS $VCF | cut -f 11 | cut -f 1 -d :`
else
GT_SANGUE=`grep -w -e ^$CHR$'\t'$POS $VCF | cut -f 10 | cut -f 1 -d :`
GT_TECIDO=`grep -w -e ^$CHR$'\t'$POS $VCF | cut -f 11 | cut -f 1 -d :`
fi

cobertura_TECIDO=`grep -w -e ^$CHR$'\t'$POS $VCF | cut -f 11 | cut -f 3 -d :`
cobertura_SANGUE=`grep -w -e ^$CHR$'\t'$POS $VCF | cut -f 10 | cut -f 3 -d :`

alteracao=`grep -w -e ^$CHR$'\t'$POS $VCF | cut -f 4`

coberturaALT_TECIDO=`grep -w -e ^$CHR$'\t'$POS $VCF | cut -f 11 | cut -f 5 -d: | cut -f 1 -d ,`
coberturaALT_SANGUE=`grep -w -e ^$CHR$'\t'$POS $VCF | cut -f 10 | cut -f 5 -d: | cut -f 1 -d ,`

if [ "$alteracao" == "A" ]
then
coberturaALT_TECIDO=`grep -w -e ^$CHR$'\t'$POS $VCF | cut -f 11 | cut -f 5 -d: | cut -f 1 -d ,`
coberturaALT_SANGUE=`grep -w -e ^$CHR$'\t'$POS $VCF | cut -f 10 | cut -f 5 -d: | cut -f 1 -d ,`
elif [ "$alteracao" == "C" ]
then
coberturaALT_TECIDO=`grep -w -e ^$CHR$'\t'$POS $VCF | cut -f 11 | cut -f 5 -d: | cut -f 2 -d ,`
coberturaALT_SANGUE=`grep -w -e ^$CHR$'\t'$POS $VCF | cut -f 10 | cut -f 5 -d: | cut -f 2 -d ,`
elif [ "$alteracao" == "G" ]
then
coberturaALT_TECIDO=`grep -w -e ^$CHR$'\t'$POS $VCF | cut -f 11 | cut -f 5 -d: | cut -f 3 -d ,`
coberturaALT_SANGUE=`grep -w -e ^$CHR$'\t'$POS $VCF | cut -f 10 | cut -f 5 -d: | cut -f 3 -d ,`
elif [ "$alteracao" == "T" ]
then
coberturaALT_TECIDO=`grep -w -e ^$CHR$'\t'$POS $VCF | cut -f 11 | cut -f 5 -d: | cut -f 4 -d ,`
coberturaALT_SANGUE=`grep -w -e ^$CHR$'\t'$POS $VCF | cut -f 10 | cut -f 5 -d: | cut -f 4 -d ,`
fi

MOSAICO_TECIDO=`bc <<< "100 - ($coberturaALT_TECIDO * 100 / ($cobertura_TECIDO))"`
MOSAICO_SANGUE=`bc <<< "100 - ($coberturaALT_SANGUE * 100 / ($cobertura_SANGUE))"`

variant_status_SANGUE=`grep -w -e ^$CHR$'\t'$POS $VCF | cut -f 10 | cut -f 12 -d:`
variant_status_TECIDO=`grep -w -e ^$CHR$'\t'$POS $VCF | cut -f 11 | cut -f 12 -d:`

echo $otherinfo $IMPACT $SYMBOL $BIOTYPE $GMAF $SIFT $PolyPhen $GT_TECIDO $cobertura_TECIDO $MOSAICO_TECIDO $variant_status_TECIDO $GT_SANGUE $cobertura_SANGUE $MOSAICO_SANGUE $variant_status_SANGUE
done < $VCF.CABECALHO
rm $VCF.CABECALHO
}

LC_ALL=C locale

for vcf in `ls *MTOR.vcf`
do
	# VEP $vcf
	echo PROCESSING $vcf ... &&
	VEPtoTable $vcf > $vcf.txt 2> $vcf.txt.log &&
	cat ~murilo/Dropbox/sh_scripts/VEPtoTable_SomaticSniper.header $vcf.txt > $vcf.txt2 &&
	mv  $vcf.txt2 $vcf.txt &&
	echo CONCLUDED WITH $vcf 
done




echo $otherinfo $IMPACT $SYMBOL $BIOTYPE $GMAF $SIFT $PolyPhen $GT_TECIDO $cobertura_TECIDO $MOSAICO_TECIDO $variant_status_TECIDO $GT_SANGUE $cobertura_SANGUE $MOSAICO_SANGUE $variant_status_SANGUE

grep -w -e ^$CHR$'\t'$POS $VCF

samtools depth ../../../../exoma23_res/bams/Sample_977_14.merged.bam.realn.bam ../../../../exoma23_res/bams/Sample_g120.merged.bam.realn.bam -r $CHR:`bc <<< "$POS-1"`-`bc <<< "$POS +1"`


