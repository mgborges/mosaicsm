TECIDO=/iscsi/murilo/exoma_4_23_25/bams/52313.merged.bam.realn.bam
SANGUE=/iscsi/murilo/exoma_4_23_25/bams/G155.merged.bam.realn.bam

pasta1=`basename $SANGUE | cut -f 1 -d _`
pasta2=`basename $TECIDO | cut -f 1 -d _`
PASTA=`echo $pasta1-$pasta2`

VCF=/home/bioinf/exoma/exoma_4_23_25/results/G155.52313.vcf
BED=/home/murilo/Dropbox/Nextera/Agillent/SureSelectXTAllexon_v6_HG38.bed 
THREADS=5

REF="/home/bioinf/ref/Homo_sapiens/NCBI/GRCh38Decoy/Sequence/BWAIndex/genome.fa" 

# MosaicHunter
java -jar /home/bioinf/exoma/testMosaicismTools/bin/MosaicHunter/build/mosaichunter.jar exome -P input_file=$TECIDO -P reference_file=$REF -P output_dir=/home/bioinf/exoma/testMosaicismTools/results/MosaicHunter/$PASTA -P depth_filter.min_depth=200 -P mode=paired_naive -P control_bam_file=$SANGUE

# LocHap
# NÃO ACEITA GZIP!!!
mkdir /home/bioinf/exoma/testMosaicismTools/results/LocHap/
samtools index $TECIDO
samtools index $SANGUE
samplename=`echo $pasta2 | cut -f 1 -d .`
/home/bioinf/exoma/testMosaicismTools/bin/LocHap-release-v2.0/bin/LocHap --vcf $VCF --bam $TECIDO --sample $samplename --out /home/bioinf/exoma/testMosaicismTools/results/LocHap/$pasta2 

# LUMPY
#mkdir /home/bioinf/exoma/testMosaicismTools/results/lumpy
#lumpyexpress -B $TECIDO,$SANGUE -o /home/bioinf/exoma/testMosaicismTools/results/lumpy/$PASTA.vcf 

# JointSNVMix
cd /home/bioinf/exoma/testMosaicismTools/bin/JointSNVMix 
mkdir /home/bioinf/exoma/testMosaicismTools/results/JointSNVMix/
jsm.py classify --model snvmix2 --out_file /home/bioinf/exoma/testMosaicismTools/results/JointSNVMix/$PASTA.txt $REF $SANGUE $TECIDO 

# VARSCAN2
mkdir /home/bioinf/exoma/testMosaicismTools/results/varscan/ 
cd /home/bioinf/exoma/testMosaicismTools/results/varscan/ 

samtools mpileup -f $REF $SANGUE $TECIDO > $PASTA.pileup 

java -jar /home/bioinf/exoma/testMosaicismTools/bin/varscan/VarScan.v2.4.3.jar somatic  $PASTA.pileup $PASTA.somatic --mpileup 1 &&
java -jar /home/bioinf/exoma/testMosaicismTools/bin/varscan/VarScan.v2.4.3.jar copynumber $PASTA.pileup $PASTA.copynumber --mpileup 1 
rm $PASTA.pileup 

# ADTEX
#mkdir /home/bioinf/exoma/testMosaicismTools/results/ADTEx 
#cd /home/bioinf/exoma/testMosaicismTools/results/ADTEx 
#python /home/bioinf/exoma/testMosaicismTools/bin/ADTEx.v.2.0/ADTEx.py --normal $SANGUE --tumor $TECIDO -o $PASTA --minReadDepth 10 -b $BED

# MUTECT2
mkdir /home/bioinf/exoma/testMosaicismTools/results/MUTECT2/
java -jar /opt/GenomeAnalysisTK.jar -T MuTect2 -R $REF -I:tumor $TECIDO -I:normal $SANGUE -o /home/bioinf/exoma/testMosaicismTools/results/MUTECT2/$PASTA.vcf

# somatic-sniper
mkdir /home/bioinf/exoma/testMosaicismTools/results/SomaticSniper
/home/bioinf/exoma/testMosaicismTools/bin/somatic-sniper/build/bin/bam-somaticsniper -F vcf -Q 40 -L -G -f $REF $TECIDO $SANGUE /home/bioinf/exoma/testMosaicismTools/results/SomaticSniper/$PASTA.vcf 

# STRELKA
STRELKA_INSTALL_DIR=/home/bioinf/exoma/testMosaicismTools/bin/strelka-2.9.2.centos6_x86_64
${STRELKA_INSTALL_DIR}/bin/configureStrelkaSomaticWorkflow.py \
    --normalBam $SANGUE \
    --tumorBam $TECIDO \
    --ref $REF \
    --exome \
    --runDir /home/bioinf/exoma/testMosaicismTools/results/STRELKA/$PASTA.somatic
# execution on a single local machine with 20 parallel jobs
/home/bioinf/exoma/testMosaicismTools/results/STRELKA/$PASTA.somatic/runWorkflow.py -m local -j $THREADS

${STRELKA_INSTALL_DIR}/bin/configureStrelkaGermlineWorkflow.py \
    --bam $SANGUE \
    --bam $TECIDO \
    --ref $REF \
    --exome \
    --runDir /home/bioinf/exoma/testMosaicismTools/results/STRELKA/$PASTA.germline
# execution on a single local machine with 20 parallel jobs
/home/bioinf/exoma/testMosaicismTools/results/STRELKA/$PASTA.germline/runWorkflow.py -m local -j $THREADS


# EBCALL
mkdir /home/bioinf/exoma/testMosaicismTools/results/EBCALL
cd /home/bioinf/exoma/testMosaicismTools/bin/EBCall
./ebCall_v2.sh $TECIDO $SANGUE /home/bioinf/exoma/testMosaicismTools/results/EBCALL/$PASTA /home/bioinf/exoma/testMosaicismTools/bin/EBCall/painelDeNormais.txt

