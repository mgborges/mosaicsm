REF="/home/bioinf/ref/Homo_sapiens/NCBI/GRCh38Decoy/Sequence/BWAIndex/genome.fa" 
THREADS=1
cd /home/bioinf/paineis/MTOR_set_18/mosaicos

function STRELKA
{
	SANGUE=$1
	TECIDO=$2
	pasta1=`basename $SANGUE | cut -f 1 -d _`
	pasta2=`basename $TECIDO | cut -f 1 -d _`
	PASTA=`echo $pasta2-$pasta1`

	STRELKA_INSTALL_DIR=/home/bioinf/exoma/testMosaicismTools/bin/strelka-2.9.2.centos6_x86_64
	${STRELKA_INSTALL_DIR}/bin/configureStrelkaSomaticWorkflow.py \
	    --normalBam $SANGUE \
	    --tumorBam $TECIDO \
	    --ref $REF \
	    --exome \
	    --runDir /home/bioinf/paineis/MTOR_set_18/mosaicos/STRELKA/$PASTA.somatic
	# execution on a single local machine with 20 parallel jobs
	/home/bioinf/paineis/MTOR_set_18/mosaicos/STRELKA/$PASTA.somatic/runWorkflow.py -m local -j $THREADS


vcf-concat /home/bioinf/paineis/MTOR_set_18/mosaicos/STRELKA/$PASTA.somatic/results/variants/somatic.indels.vcf.gz /home/bioinf/paineis/MTOR_set_18/mosaicos/STRELKA/$PASTA.somatic/results/variants/somatic.snvs.vcf.gz | bcftools view -f PASS | java -jar /opt/snpEff/SnpSift.jar filter "(DP>200)" > /home/bioinf/paineis/MTOR_set_18/mosaicos/STRELKA/$PASTA.filtered.vcf

}

#mkdir STRELKA
cd STRELKA

STRELKA ../../bams/74212.MERGED.bam.realn.bam ../../bams/P20.MERGED.bam.realn.bam &
STRELKA ../../bams/54112.MERGED.bam.realn.bam ../../bams/P19.MERGED.bam.realn.bam &
STRELKA ../../bams/52313.MERGED.bam.realn.bam ../../bams/G155.MERGED.bam.realn.bam &
STRELKA ../../bams/116015.MERGED.bam.realn.bam ../../bams/G150.MERGED.bam.realn.bam &
STRELKA ../../bams/39015.MERGED.bam.realn.bam ../../bams/G133.MERGED.bam.realn.bam &
STRELKA ../../bams/89114.MERGED.bam.realn.bam ../../bams/G118.MERGED.bam.realn.bam &
STRELKA ../../bams/81410.MERGED.bam.realn.bam ../../bams/P02.MERGED.bam.realn.bam &
STRELKA ../../bams/17715.MERGED.bam.realn.bam ../../bams/G129.MERGED.bam.realn.bam &
STRELKA ../../bams/97714.MERGED.bam.realn.bam ../../bams/G120.MERGED.bam.realn.bam &
STRELKA ../../bams/21315.MERGED.bam.realn.bam ../../bams/G9.MERGED.bam.realn.bam &
STRELKA ../../bams/G157.MERGED.bam.realn.bam ../../bams/46916.MERGED.bam.realn.bam &
STRELKA ../../bams/15516.MERGED.bam.realn.bam ../../bams/G125.MERGED.bam.realn.bam &



function VEP_VVP
{
	VCF=$1
	/home/murilo/VVP-pub/vt/vt decompose -s $VCF -o $VCF.decomp.vcf
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

	java -jar /home/murilo/snpEff/SnpSift.jar extractFields $VCF CHROM POS ALT REF "GEN[*].GT" "GEN[*].DP" "FILTER" "GEN[*]" -e "." -s ";" > $VCF.tmp

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
	GT_INFO=`grep -w -e ^$CHR$'\t'$POS $VCF.tmp | cut -f 8`
	elif [ "$verificator" != "0" ] && [ "$ALL" == "-" ]
	then
	POS=`bc <<< "$POS - 1"`
	GT=`grep -w -e ^$CHR$'\t'$POS $VCF.tmp | cut -f 5`
	DP=`grep -w -e ^$CHR$'\t'$POS $VCF.tmp | cut -f 6`
	GT_INFO=`grep -w -e ^$CHR$'\t'$POS $VCF.tmp | cut -f 8`
	elif [ "$verificator" == "0" ] && [ "$ALL" == "-" ]
	then
	POS=`bc <<< "$POS - 1"`
	GT=`grep -w -e ^$CHR$'\t'$POS $VCF.tmp | cut -f 5`
	DP=`grep -w -e ^$CHR$'\t'$POS $VCF.tmp | cut -f 6`
	GT_INFO=`grep -w -e ^$CHR$'\t'$POS $VCF.tmp | cut -f 8`
	else
	GT=`grep -w -e ^$CHR$'\t'$POS $VCF.tmp | cut -f 5`
	DP=`grep -w -e ^$CHR$'\t'$POS $VCF.tmp | cut -f 6`
	GT_INFO=`grep -w -e ^$CHR$'\t'$POS $VCF.tmp | cut -f 8`
	fi

ACGT_SANGUE=`echo $GT_INFO | cut -f 1 -d ';' | cut -f 5-8 -d : | sed 's/:/,/g' | cut -f 1,3,5,7 -d ,`
ACGT_TECIDO=`echo $GT_INFO | cut -f 2 -d ';' | cut -f 5-8 -d : | sed 's/:/,/g' | cut -f 1,3,5,7 -d ,`

	cobertura_TECIDO=`echo $DP | cut -f 2 -d ';'`
	cobertura_SANGUE=`echo $DP | cut -f 1 -d ';'`

	echo $otherinfo $IMPACT $SYMBOL $BIOTYPE $GMAF $SIFT $PolyPhen $FILTER $cobertura_TECIDO $ACGT_TECIDO $cobertura_SANGUE $ACGT_SANGUE

	done < $VCF.CABECALHO
	rm $VCF.tmp $VCF.CABECALHO
}

function EXECUTE_VEP_VVP
{
	vcf=$1
	PATHWAY=$2
	VEP_VVP $vcf && 
	VEPtoTable $vcf.decomp.VEP.vcf | grep -e shift -e missense_variant -e nonsense -e splice -e stop > $vcf.$PATHWAY.tmp && 
	cat ~murilo/Dropbox/sh_scripts/VEPtoTable_VVP_STRELKA.header $vcf.$PATHWAY.tmp > $vcf.$PATHWAY.filtered.xls && 
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

for VCF_FILES in `ls *vcf`
do
/home/murilo/ensembl-vep/vep -i $VCF_FILES -o $VCF.VEP.vcf --cache  --pubmed --sift b --polyphen b --ccds --uniprot --symbol --numbers --domains --regulatory --canonical --protein --biotype --gene_phenotype --af --variant_class --force_overwrite --vcf --merged
done
