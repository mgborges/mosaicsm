REF="/home/bioinf/ref/Homo_sapiens/NCBI/GRCh38Decoy/Sequence/BWAIndex/genome.fa" 
THREADS=1
cd /home/bioinf/paineis/MTOR_set_18/mosaicos

function JOINTSNV
{
	SANGUE=$1
	TECIDO=$2
	pasta1=`basename $SANGUE | cut -f 1 -d _`
	pasta2=`basename $TECIDO | cut -f 1 -d _`
	PASTA=`echo $pasta2-$pasta1`

cd /home/bioinf/exoma/testMosaicismTools/bin/JointSNVMix 
echo indexing $SANGUE
#samtools index $SANGUE
echo indexing $TECIDO
#samtools index $TECIDO
#jsm.py classify --model snvmix2 --out_file /home/bioinf/paineis/MTOR_set_18/mosaicos/JOINTSNV/$PASTA.txt $REF $SANGUE $TECIDO 

cat /home/bioinf/paineis/MTOR_set_18/mosaicos/JOINTSNV/$PASTA.txt | awk '{if(($5+$6)>200)print}' | awk '{if(($7+$8)>200)print}' | awk '{if(($8*100/($7+$8))<20)print}' | awk '{if($8>2)print}' | awk '{if(($6*100/($5+$6))<1)print}' | awk '{if($9>0.80)print}' | awk '{if($9<1)print}' > /home/bioinf/paineis/MTOR_set_18/mosaicos/JOINTSNV/$PASTA.filtered.tmp

less /home/bioinf/paineis/MTOR_set_18/mosaicos/JOINTSNV/$PASTA.filtered.tmp | awk '{print $1 "\t" $2 "\t" $2 "\t" $3 "/"$4 "\t1"}' > /home/bioinf/paineis/MTOR_set_18/mosaicos/JOINTSNV/$PASTA.filtered.TOVEP.tmp

}

#mkdir JOINTSNV
cd JOINTSNV

JOINTSNV /home/bioinf/paineis/MTOR_set_18/bams/74212.MERGED.bam.realn.bam /home/bioinf/paineis/MTOR_set_18/bams/P20.MERGED.bam.realn.bam &
JOINTSNV /home/bioinf/paineis/MTOR_set_18/bams/54112.MERGED.bam.realn.bam /home/bioinf/paineis/MTOR_set_18/bams/P19.MERGED.bam.realn.bam &
JOINTSNV /home/bioinf/paineis/MTOR_set_18/bams/52313.MERGED.bam.realn.bam /home/bioinf/paineis/MTOR_set_18/bams/G155.MERGED.bam.realn.bam &
JOINTSNV /home/bioinf/paineis/MTOR_set_18/bams/116015.MERGED.bam.realn.bam /home/bioinf/paineis/MTOR_set_18/bams/G150.MERGED.bam.realn.bam &
JOINTSNV /home/bioinf/paineis/MTOR_set_18/bams/39015.MERGED.bam.realn.bam /home/bioinf/paineis/MTOR_set_18/bams/G133.MERGED.bam.realn.bam &
JOINTSNV /home/bioinf/paineis/MTOR_set_18/bams/89114.MERGED.bam.realn.bam /home/bioinf/paineis/MTOR_set_18/bams/G118.MERGED.bam.realn.bam &
JOINTSNV /home/bioinf/paineis/MTOR_set_18/bams/81410.MERGED.bam.realn.bam /home/bioinf/paineis/MTOR_set_18/bams/P02.MERGED.bam.realn.bam &
JOINTSNV /home/bioinf/paineis/MTOR_set_18/bams/17715.MERGED.bam.realn.bam /home/bioinf/paineis/MTOR_set_18/bams/G129.MERGED.bam.realn.bam &
JOINTSNV /home/bioinf/paineis/MTOR_set_18/bams/97714.MERGED.bam.realn.bam /home/bioinf/paineis/MTOR_set_18/bams/G120.MERGED.bam.realn.bam &
JOINTSNV /home/bioinf/paineis/MTOR_set_18/bams/21315.MERGED.bam.realn.bam /home/bioinf/paineis/MTOR_set_18/bams/G9.MERGED.bam.realn.bam &
JOINTSNV /home/bioinf/paineis/MTOR_set_18/bams/G157.MERGED.bam.realn.bam /home/bioinf/paineis/MTOR_set_18/bams/46916.MERGED.bam.realn.bam &
JOINTSNV /home/bioinf/paineis/MTOR_set_18/bams/15516.MERGED.bam.realn.bam /home/bioinf/paineis/MTOR_set_18/bams/G125.MERGED.bam.realn.bam &



function VEP_VVP
{
	VCF=$1
	
	/home/murilo/ensembl-vep/vep -i $VCF -o $VCF.VEP.vcf --cache  --pubmed --sift b --polyphen b --ccds --uniprot --symbol --numbers --domains --regulatory --canonical --protein --biotype --gene_phenotype --af --variant_class --force_overwrite --vcf --merged

	/home/murilo/ensembl-vep/vep -i $VCF -o $VCF.VEP --cache  --pubmed --sift b --polyphen b --ccds --uniprot --symbol --numbers --domains --regulatory --canonical --protein --biotype --gene_phenotype --af --variant_class --force_overwrite  --merged
}

function VEPtoTable
{
	VCF=$1
	VEP=`echo $VCF | sed 's/.vcf$//g'`
	grep -v ^# $VEP > $VCF.CABECALHO
	TXT=`echo $VCF | cut -f 1 -d .`


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

SSINFO=`grep -w -e ^$CHR$'\t'$POS $TXT*.txt | cut -f 5-17 | sed 's/\t/ /g'`

	echo $otherinfo $IMPACT $SYMBOL $BIOTYPE $GMAF $SIFT $PolyPhen $SSINFO

	done < $VCF.CABECALHO
	rm $VCF.CABECALHO
}

function EXECUTE_VEP_VVP
{
	vcf=$1
	PATHWAY=$2
	VEP_VVP $vcf && 
	VEPtoTable $vcf.VEP.vcf | grep -e shift -e missense_variant -e nonsense -e splice -e stop > $vcf.$PATHWAY.tmp && 
	cat ~murilo/Dropbox/sh_scripts/VEPtoTable_VVP_JOINTSNV.header $vcf.$PATHWAY.tmp > $vcf.$PATHWAY.filtered.xls && 
	rm $vcf*VEP $vcf $vcf*VEP_summary.html
	echo DONE FOR $vcf
}

# Quanto aos filtros, pode aplicar para as consequencias: frame-shift, missense, nonsense, splicing site, stop-codon

####### WHOLE EXOME

PATHWAY="ALL"

for vcf in `ls *TOVEP*.tmp`
do
EXECUTE_VEP_VVP $vcf $PATHWAY 
done
