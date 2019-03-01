REF="/home/bioinf/ref/Homo_sapiens/NCBI/GRCh38Decoy/Sequence/BWAIndex/genome.fa" 
THREADS=1
cd /home/bioinf/paineis/MTOR_set_18/mosaicos

function LOCHAP
{
	TECIDO=$1
	pasta2=`basename $TECIDO | cut -f 1 -d _`
	VCF=$2

# LocHap
# NÃO ACEITA GZIP!!!

samplename=`echo $pasta2 | cut -f 1 -d .`
/home/bioinf/exoma/testMosaicismTools/bin/LocHap-release-v2.0/bin/LocHap --vcf $VCF --bam $TECIDO --sample $samplename --out $pasta2 --igv

grep -v ^# $pasta2.igv > $pasta2.igv.tmp
grep ^# $pasta2.igv > $pasta2.igv.HEADER

while read -r line
do
DP=`echo $line | grep nTot| awk -F 'nTot=' '{print $2}' | cut -f 1 -d '<'`
if [ "$DP" -gt 200 ]
then
echo $line | sed 's/ /\t/g'
fi
done < $pasta2.igv.tmp > $pasta2.igv.lines

cat $pasta2.igv.HEADER $pasta2.igv.lines > $pasta2.igv

rm $pasta2.igv.tmp $pasta2.igv.HEADER $pasta2.igv.lines
}

#mkdir LOCHAP
cd LOCHAP

LOCHAP ../../bams/P20.MERGED.bam.realn.bam ../../results/74212-P20.vcf
LOCHAP ../../bams/P19.MERGED.bam.realn.bam ../../results/54112-P19.vcf
LOCHAP ../../bams/G155.MERGED.bam.realn.bam ../../results/52313-G155.vcf
LOCHAP ../../bams/G150.MERGED.bam.realn.bam ../../results/116015-G150.vcf
LOCHAP ../../bams/G133.MERGED.bam.realn.bam ../../results/39015-G133.vcf
LOCHAP ../../bams/G118.MERGED.bam.realn.bam ../../results/89114-G118.vcf
LOCHAP ../../bams/P02.MERGED.bam.realn.bam ../../results/81410-P02.vcf
LOCHAP ../../bams/G129.MERGED.bam.realn.bam ../../results/17715-G129.vcf
LOCHAP ../../bams/G120.MERGED.bam.realn.bam ../../results/97714-G120.vcf
LOCHAP ../../bams/G9.MERGED.bam.realn.bam ../../results/21315-G9.vcf
LOCHAP ../../bams/G157.MERGED.bam.realn.bam ../../results/46916-G157.vcf
LOCHAP ../../bams/G125.MERGED.bam.realn.bam ../../results/15516-G125.vcf

