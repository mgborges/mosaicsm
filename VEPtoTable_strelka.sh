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
GMAF=`echo $line | awk -F "AF=" '{print $2}' | awk -F ";" '{print $1}'`
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
GT=`grep -w -e ^$CHR$'\t'$POS $VCF | cut -f 5`
DP=`grep -w -e ^$CHR$'\t'$POS $VCF | cut -f 6`
elif [ "$verificator" != "0" ] && [ "$ALL" == "-" ]
then
POS=`bc <<< "$POS - 1"`
GT=`grep -w -e ^$CHR$'\t'$POS $VCF | cut -f 5`
DP=`grep -w -e ^$CHR$'\t'$POS $VCF | cut -f 6`
elif [ "$verificator" == "0" ] && [ "$ALL" == "-" ]
then
POS=`bc <<< "$POS - 1"`
GT=`grep -w -e ^$CHR$'\t'$POS $VCF | cut -f 5`
DP=`grep -w -e ^$CHR$'\t'$POS $VCF | cut -f 6`
else
GT=`grep -w -e ^$CHR$'\t'$POS $VCF | cut -f 5`
DP=`grep -w -e ^$CHR$'\t'$POS $VCF | cut -f 6`
fi

cobertura_TECIDO=`grep -w -e ^$CHR$'\t'$POS $VCF | cut -f 11 | cut -f 1 -d:`
cobertura_SANGUE=`grep -w -e ^$CHR$'\t'$POS $VCF | cut -f 10 | cut -f 1 -d:`

alteracao=`grep -w -e ^$CHR$'\t'$POS $VCF | cut -f 4`

coberturaALT_TECIDO=`grep -w -e ^$CHR$'\t'$POS $VCF | cut -f 11 | cut -f 5 -d: | cut -f 1 -d ,`
coberturaALT_SANGUE=`grep -w -e ^$CHR$'\t'$POS $VCF | cut -f 10 | cut -f 5 -d: | cut -f 1 -d ,`

MOSAICO_TECIDO=`bc <<< "100 - ($coberturaALT_TECIDO * 100 / ($cobertura_TECIDO))"`
MOSAICO_SANGUE=`bc <<< "100 - ($coberturaALT_SANGUE * 100 / ($cobertura_SANGUE))"`

echo $otherinfo $IMPACT $SYMBOL $BIOTYPE $GMAF $SIFT $PolyPhen $cobertura_TECIDO $MOSAICO_TECIDO $cobertura_SANGUE $MOSAICO_SANGUE
done < $VCF.CABECALHO
rm $VCF.CABECALHO
}

LC_ALL=C locale

###### FOR STRELKA
for vcf in `ls */results/*pass*.MTOR.vcf`
do
	# VEP $vcf
	echo PROCESSING $vcf ... &&
	VEPtoTable $vcf > $vcf.txt 2> $vcf.txt.log &&
	cat ~murilo/Dropbox/sh_scripts/VEPtoTable_strelka.header $vcf.txt > $vcf.txt2 &&
	mv  $vcf.txt2 $vcf.txt &&
	echo CONCLUDED WITH $vcf
done


