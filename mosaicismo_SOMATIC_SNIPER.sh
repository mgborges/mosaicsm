REF="/home/bioinf/ref/Homo_sapiens/NCBI/GRCh38Decoy/Sequence/BWAIndex/genome.fa" 
THREADS=1
cd /home/bioinf/paineis/MTOR_set_18/mosaicos

function SOMSNIPER
{
	SANGUE=$1
	TECIDO=$2
	pasta1=`basename $SANGUE | cut -f 1 -d _`
	pasta2=`basename $TECIDO | cut -f 1 -d _`
	PASTA=`echo $pasta2-$pasta1`

/home/bioinf/exoma/testMosaicismTools/bin/somatic-sniper/build/bin/bam-somaticsniper -F vcf -Q 25 -f $REF $TECIDO $SANGUE $PASTA.vcf -n $pasta1 -t $pasta2

java -jar /opt/snpEff/SnpSift.jar filter "(DP>200)" $PASTA.vcf > $PASTA.filtered.vcf
}

#mkdir SOMSNIPER
cd SOMSNIPER

SOMSNIPER ../../bams/74212.MERGED.bam.realn.bam ../../bams/P20.MERGED.bam.realn.bam &
SOMSNIPER ../../bams/54112.MERGED.bam.realn.bam ../../bams/P19.MERGED.bam.realn.bam &
SOMSNIPER ../../bams/52313.MERGED.bam.realn.bam ../../bams/G155.MERGED.bam.realn.bam &
SOMSNIPER ../../bams/116015.MERGED.bam.realn.bam ../../bams/G150.MERGED.bam.realn.bam &
SOMSNIPER ../../bams/39015.MERGED.bam.realn.bam ../../bams/G133.MERGED.bam.realn.bam &
SOMSNIPER ../../bams/89114.MERGED.bam.realn.bam ../../bams/G118.MERGED.bam.realn.bam &
SOMSNIPER ../../bams/81410.MERGED.bam.realn.bam ../../bams/P02.MERGED.bam.realn.bam &
SOMSNIPER ../../bams/17715.MERGED.bam.realn.bam ../../bams/G129.MERGED.bam.realn.bam &
SOMSNIPER ../../bams/97714.MERGED.bam.realn.bam ../../bams/G120.MERGED.bam.realn.bam &
SOMSNIPER ../../bams/21315.MERGED.bam.realn.bam ../../bams/G9.MERGED.bam.realn.bam &
SOMSNIPER ../../bams/G157.MERGED.bam.realn.bam ../../bams/46916.MERGED.bam.realn.bam &
SOMSNIPER ../../bams/15516.MERGED.bam.realn.bam ../../bams/G125.MERGED.bam.realn.bam &


