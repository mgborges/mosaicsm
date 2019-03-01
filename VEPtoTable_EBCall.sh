#Location Allele Gene Feature Feature_type Consequence cDNA_position CDS_position Protein_position Amino_acids Codons Existing_variation IMPACT SYMBOL BIOTYPE GMAF SIFT PolyPhen GT DP MOSAICO_TECIDO cobertura_TECIDO MOSAICO_SANGUE cobertura_SANGUE

cd /home/bioinf/exoma/exoma_25_res/results/mosaicos/EBCall/EBCall

function VEPtoTable
{
export VEP="/home/murilo/ensembl-tools-release-84/scripts/variant_effect_predictor/variant_effect_predictor.pl"

INPUT="$1/output.txt"

cp $INPUT TMP.tmp

less TMP.tmp | grep -v Chr | cut -f 1-5 > tmp.toVEP

perl $VEP -i tmp.toVEP -o tmp.VEP --cache  --pubmed --sift b --polyphen b --ccds --uniprot --symbol --numbers --domains --regulatory --canonical --protein --biotype --gene_phenotype --gmaf --variant_class --port 3337 --force_overwrite --merged --fork 4

grep -v ^# tmp.VEP > TMP.CABECALHO

rm tmp.toVEP *html *warnings* tmp.semCabecalho #tmp.VEP

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
otherinfo=`echo $line | cut -f 4-13 -d " "`

CHR=`echo $line | cut -f 2 -d " " | cut -f 1 -d ":"`
CHR="chr$CHR"

POS=`echo $line | cut -f 2 -d " " | cut -f 2 -d ":" | cut -f 1 -d "-"`

REF=`grep -w -e $CHR$'\t'$POS TMP.tmp | cut -f 4`
ALT=`grep -w -e $CHR$'\t'$POS TMP.tmp | cut -f 5`

MOSAICO_TECIDO=`grep -w -e ^$CHR$'\t'$POS TMP.tmp | cut -f 9`
MOSAICO_SANGUE=`grep -w -e ^$CHR$'\t'$POS TMP.tmp | cut -f 13`

cobertura_TECIDO=`grep -w -e ^$CHR$'\t'$POS TMP.tmp | cut -f 10`
cobertura_SANGUE=`grep -w -e ^$CHR$'\t'$POS TMP.tmp | cut -f 14`

echo $CHR $POS $REF $ALT $otherinfo $IMPACT $SYMBOL $BIOTYPE $GMAF $SIFT $PolyPhen $MOSAICO_TECIDO $cobertura_TECIDO $MOSAICO_SANGUE $cobertura_SANGUE >> tmp.semCabecalho
done < TMP.CABECALHO

cat /home/murilo/Dropbox/sh_scripts/VEPtoTable_EBCall.header tmp.semCabecalho > $folder.EBCall.txt
rm TMP.CABECALHO tmp.VEP TMP.tmp tmp.semCabecalho
}

for folder in `ls`
do
	echo PROCESSING $folder ...
	VEPtoTable $folder
done
