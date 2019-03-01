REF="/home/bioinf/ref/Homo_sapiens/NCBI/GRCh38Decoy/Sequence/BWAIndex/genome.fa" 
cd /home/bioinf/paineis/MTOR_set_18/mosaicos

function callMUTECT2
{
	SANGUE=$1
	TECIDO=$2
	pasta1=`basename $SANGUE | cut -f 1 -d _`
	pasta2=`basename $TECIDO | cut -f 1 -d _`
	PASTA=`echo $pasta2-$pasta1`
	java -jar /opt/GenomeAnalysisTK.jar -T MuTect2 -R $REF -I:tumor $TECIDO -I:normal $SANGUE -o $PASTA.vcf
}
#mkdir MUTECT2
cd MUTECT2

callMUTECT2 /home/bioinf/paineis/MTOR_set_18/bams/74212.MERGED.bam.realn.bam /home/bioinf/paineis/MTOR_set_18/bams/P20.MERGED.bam.realn.bam &
callMUTECT2 /home/bioinf/paineis/MTOR_set_18/bams/54112.MERGED.bam.realn.bam /home/bioinf/paineis/MTOR_set_18/bams/P19.MERGED.bam.realn.bam &
callMUTECT2 /home/bioinf/paineis/MTOR_set_18/bams/52313.MERGED.bam.realn.bam /home/bioinf/paineis/MTOR_set_18/bams/G155.MERGED.bam.realn.bam &
callMUTECT2 /home/bioinf/paineis/MTOR_set_18/bams/116015.MERGED.bam.realn.bam /home/bioinf/paineis/MTOR_set_18/bams/G150.MERGED.bam.realn.bam &
callMUTECT2 /home/bioinf/paineis/MTOR_set_18/bams/39015.MERGED.bam.realn.bam /home/bioinf/paineis/MTOR_set_18/bams/G133.MERGED.bam.realn.bam &
callMUTECT2 /home/bioinf/paineis/MTOR_set_18/bams/89114.MERGED.bam.realn.bam /home/bioinf/paineis/MTOR_set_18/bams/G118.MERGED.bam.realn.bam &
callMUTECT2 /home/bioinf/paineis/MTOR_set_18/bams/81410.MERGED.bam.realn.bam /home/bioinf/paineis/MTOR_set_18/bams/P02.MERGED.bam.realn.bam &
callMUTECT2 /home/bioinf/paineis/MTOR_set_18/bams/17715.MERGED.bam.realn.bam /home/bioinf/paineis/MTOR_set_18/bams/G129.MERGED.bam.realn.bam &
callMUTECT2 /home/bioinf/paineis/MTOR_set_18/bams/97714.MERGED.bam.realn.bam /home/bioinf/paineis/MTOR_set_18/bams/G120.MERGED.bam.realn.bam &
callMUTECT2 /home/bioinf/paineis/MTOR_set_18/bams/21315.MERGED.bam.realn.bam /home/bioinf/paineis/MTOR_set_18/bams/G9.MERGED.bam.realn.bam &
callMUTECT2 /home/bioinf/paineis/MTOR_set_18/bams/G157.MERGED.bam.realn.bam /home/bioinf/paineis/MTOR_set_18/bams/46916.MERGED.bam.realn.bam &
callMUTECT2 /home/bioinf/paineis/MTOR_set_18/bams/15516.MERGED.bam.realn.bam /home/bioinf/paineis/MTOR_set_18/bams/G125.MERGED.bam.realn.bam &


function VEP_VVP
{
	VCF=$1
	/home/murilo/VVP-pub/vt/vt decompose -s $VCF -o $VCF.decomp.vcf
	mv $VCF.decomp.vcf $VCF.decomp.vcf_tmp
	less $VCF.decomp.vcf_tmp | grep -v  alt_allele_in_normal  | grep -v germline_risk | grep -v str_contraction > $VCF.decomp.vcf
	rm $VCF.decomp.vcf_tmp
	vcftools --vcf $VCF.decomp.vcf --recode --out $VCF
	mv $VCF.recode.vcf $VCF.decomp.vcf

	/home/murilo/ensembl-vep/vep -i $VCF.decomp.vcf -o $VCF.decomp.VEP.vcf --cache  --pubmed --sift b --polyphen b --ccds --uniprot --symbol --numbers --domains --regulatory --canonical --protein --biotype --gene_phenotype --af --variant_class --force_overwrite --vcf --merged

	/home/murilo/ensembl-vep/vep -i $VCF.decomp.vcf -o $VCF.decomp.VEP --cache  --pubmed --sift b --polyphen b --ccds --uniprot --symbol --numbers --domains --regulatory --canonical --protein --biotype --gene_phenotype --af --variant_class --force_overwrite  --merged

	cat $VCF.decomp.VEP.vcf | sed 's/^chr//' > $VCF.tmp
	rm $VCF.decomp.vcf

	/home/murilo/VVP-pub/VVP -i $VCF.tmp -d /home/murilo/VVP-pub/gnomad.062717.build -v CSQ,4,6,1,15 1> $VCF.vvp.out

	cat $VCF.vvp.out | sed 's/\./,/g' > $VCF.vvp.out.xls
	rm $VCF.tmp $VCF.vvp.out 
}

function VEPtoTable
{
	VCF=$1
	VEP=`echo $VCF | sed 's/.vcf$//g'`
	grep -v ^# $VEP > $VCF.CABECALHO
	VVP=`echo $VCF | sed 's/.decomp.VEP.vcf$/.vvp.out.xls/g'`

	java -jar /home/murilo/snpEff/SnpSift.jar extractFields $VCF CHROM POS ALT REF "GEN[*].GT" "GEN[*].AD" "FILTER" -e "." -s ";" > $VCF.tmp

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
	TRANSCRITO=`echo $otherinfo | cut -f 5 -d " "`

	CHR=`echo $line | cut -f 2 -d " " | cut -f 1 -d ":"`
	POS=`echo $line | cut -f 2 -d " " | cut -f 2 -d ":" | cut -f 1 -d "-"`
	ALL=`echo $line | cut -f 3 -d " "`

	FILTER=`grep -w -e ^$CHR$'\t'$POS $VCF | cut -f 7`

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

	REF_TECIDO=`echo $DP | cut -f 1 -d ';' | cut -f 1 -d ','`
	ALT_TECIDO=`echo $DP | cut -f 1 -d ';' | cut -f 2 -d ','`
	REF_SANGUE=`echo $DP | cut -f 2 -d ';' | cut -f 1 -d ','`
	ALT_SANGUE=`echo $DP | cut -f 2 -d ';' | cut -f 2 -d ','`

	MOSAICO_TECIDO=`bc <<< "$ALT_TECIDO * 100 / ($REF_TECIDO + $ALT_TECIDO)"`
	MOSAICO_SANGUE=`bc <<< "$ALT_SANGUE * 100 / ($REF_SANGUE + $ALT_SANGUE)"`

	cobertura_TECIDO=`bc <<< "$REF_TECIDO + $ALT_TECIDO"`
	cobertura_SANGUE=`bc <<< "$REF_SANGUE + $ALT_SANGUE"`

	VVP_SCORE=`grep $TRANSCRITO $VVP | grep $POS -m 1 | cut -f 8,13,18 | sed 's/\t/ /g'`

	echo $otherinfo $IMPACT $SYMBOL $BIOTYPE $GMAF $SIFT $PolyPhen $VVP_SCORE $FILTER $GT $DP $cobertura_TECIDO $cobertura_SANGUE $MOSAICO_TECIDO $MOSAICO_SANGUE

	done < $VCF.CABECALHO
	rm $VCF.tmp $VCF.CABECALHO
}

function EXECUTE_VEP_VVP
{
	vcf=$1
	PATHWAY=$2
	VEP_VVP $vcf && 
	VEPtoTable $vcf.decomp.VEP.vcf | grep -e shift -e missense_variant -e nonsense -e splice -e stop > $vcf.$PATHWAY.tmp && 
	cat ~murilo/Dropbox/sh_scripts/VEPtoTable_VVP_MUTECT2.header $vcf.$PATHWAY.tmp > $vcf.$PATHWAY.filtered.xls && 
	rm $vcf.$PATHWAY.tmp $vcf*decomp* $vcf*vvp.out.xls &&
	echo DONE FOR $vcf
}

# Quanto aos filtros, pode aplicar para as consequencias: frame-shift, missense, nonsense, splicing site, stop-codon

####### WHOLE EXOME

PATHWAY="ALL"

for vcf in `ls *vcf`
do
EXECUTE_VEP_VVP $vcf $PATHWAY &
done

for xlsfile in `ls *.ALL.filtered.xls`
do
output=`echo $xlsfile | sed 's/filtered/filtered_DP_MOSAIC/g'`
cat $xlsfile | awk -F" " '{if($26>200)print}' | awk -F" " '{if($27>200)print}' | awk -F" " '{if($28<20)print}' > $output.tmp
cat ~murilo/Dropbox/sh_scripts/VEPtoTable_VVP_MUTECT2.header $output.tmp > $output
rm $output.tmp $xlsfile
done

for VCF_FILES in `ls *vcf`
do
/home/murilo/ensembl-vep/vep -i $VCF_FILES -o $VCF.VEP.vcf --cache  --pubmed --sift b --polyphen b --ccds --uniprot --symbol --numbers --domains --regulatory --canonical --protein --biotype --gene_phenotype --af --variant_class --force_overwrite --vcf --merged
done
