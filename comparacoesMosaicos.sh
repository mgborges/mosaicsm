#cd /home/bioinf/exoma/exoma_25_res/results/mosaicos
cd /home/bioinf/paineis/MTOR_set_17/results/mosaicos
cd /home/bioinf/paineis/MTOR_set_17/results/mosaicos/


####################### MTOR

for file in `ls SomaticSniper/*txt`
do

SAMPLE=`echo $file | cut -f 2 -d / | cut -f 1 -d -`
NAME=`echo $file | cut -f 2 -d / | cut -f 1 -d .`

MUTECT=`ls MUTECT2/*$SAMPLE*VEP`

SOMSNIPER=`ls SomaticSniper/*$SAMPLE*VEP`

STRELKA_SNV=`ls STRELKA/*/results/all*snv*$SAMPLE*VEP`
STRELKA_INDEL=`ls STRELKA/*/results/all*indel*$SAMPLE*VEP`

#EBCALL=`ls EBCall/EBCall*_TAU/*$SAMPLE*MTOR*txt`

less $MUTECT |  grep -v ^# | cut -f 2 | sort | uniq | grep -v Uploaded_variation > $NAME.MTOR.MUTECT

less $SOMSNIPER | grep -v ^# | cut -f 2 | sort | uniq | grep -v Uploaded_variation > $NAME.MTOR.SOMSNIPER

less $STRELKA_SNV | grep -v ^# | cut -f 2 | sort | uniq | grep -v Uploaded_variation > $SAMPLE.TMP
less $STRELKA_INDEL | grep -v ^# | cut -f 2 | sort | uniq | grep -v Uploaded_variation > $SAMPLE.TMP2
cat $SAMPLE.TMP $SAMPLE.TMP2 | sort | uniq > $NAME.MTOR.STRELKA
rm $SAMPLE.TMP $SAMPLE.TMP2

#less $EBCALL | grep -v intron_variant | grep -v ^CHR | cut -f 1-4 -d " " | sed 's/chr//g' | awk -F " " '{ print $1 "_" $2 "_" $3 "/" $4}' | sort | uniq >  $NAME.NO_INTRON.MTOR.EBCALL

#less $EBCALL | grep -v ^CHR | cut -f 1-4 -d " " | sed 's/chr//g' | awk -F " " '{ print $1 "_" $2 "_" $3 "/" $4}' | sort | uniq >  $NAME.MTOR.EBCALL

done


####################### TAU

cd /home/bioinf/exoma/exoma_25_res/results/mosaicos

for file in `ls SomaticSniper/*MT*R*txt`
do

SAMPLE=`echo $file | cut -f 2 -d / | cut -f 1 -d -`
NAME=`echo $file | cut -f 2 -d / | cut -f 1 -d .`

EBCALL=`ls EBCall/EBCall_M*OR_TAU/*$SAMPLE*TAU*txt`


MUTECT=`ls mutect2/*$SAMPLE*TAU*VEP`

SOMSNIPER=`ls SomaticSniper/*$SAMPLE*TAU*VEP`

STRELKA_SNV=`ls strelka/*/all*snv*$SAMPLE*TAU*VEP`
STRELKA_INDEL=`ls strelka/*/all*indel*$SAMPLE*TAU*VEP`

EBCALL=`ls EBCall/EBCall*_TAU/*$SAMPLE*TAU*txt`

less $MUTECT | grep -v intron_variant | grep -v ^# | cut -f 1 | sort | uniq | grep -v Uploaded_variation > $NAME.NO_INTRON.TAU.MUTECT
less $MUTECT  | grep -v ^# | cut -f 1 | sort | uniq | grep -v Uploaded_variation > $NAME.TAU.MUTECT

less $SOMSNIPER | grep -v intron_variant | grep -v ^# | cut -f 1 | sort | uniq | grep -v Uploaded_variation > $NAME.NO_INTRON.TAU.SOMSNIPER
less $SOMSNIPER | grep -v ^# | cut -f 1 | sort | uniq | grep -v Uploaded_variation > $NAME.TAU.SOMSNIPER

less $STRELKA_SNV | grep -v intron_variant | grep -v ^# | cut -f 1 | sort | uniq | grep -v Uploaded_variation > $SAMPLE.TMP
less $STRELKA_INDEL | grep -v intron_variant | grep -v ^# | cut -f 1 | sort | uniq | grep -v Uploaded_variation > $SAMPLE.TMP2
cat $SAMPLE.TMP $SAMPLE.TMP2 | sort | uniq > $NAME.NO_INTRON.TAU.STRELKA
rm $SAMPLE.TMP $SAMPLE.TMP2

less $STRELKA_SNV | grep -v ^# | cut -f 1 | sort | uniq | grep -v Uploaded_variation > $SAMPLE.TMP
less $STRELKA_INDEL | grep -v ^# | cut -f 1 | sort | uniq | grep -v Uploaded_variation > $SAMPLE.TMP2
cat $SAMPLE.TMP $SAMPLE.TMP2 | sort | uniq > $NAME.TAU.STRELKA
rm $SAMPLE.TMP $SAMPLE.TMP2

less $EBCALL | grep -v intron_variant | grep -v ^CHR | cut -f 1-4 -d " " | sed 's/chr//g' | awk -F " " '{ print $1 "_" $2 "_" $3 "/" $4}' | sort | uniq >  $NAME.NO_INTRON.TAU.EBCALL

less $EBCALL | grep -v ^CHR | cut -f 1-4 -d " " | sed 's/chr//g' | awk -F " " '{ print $1 "_" $2 "_" $3 "/" $4}' | sort | uniq >  $NAME.TAU.EBCALL

done




####################### EXOMA COMPLETO

for file in `ls SOMATICSNIPER/*txt`
do

SAMPLE=`echo $file | cut -f 2 -d / | cut -f 1 -d -`
NAME=`echo $file | cut -f 2 -d / | cut -f 1 -d .`

MUTECT=`ls mutect2/*$SAMPLE*VEP | grep -v MTOR | grep -v TAU`

SOMSNIPER=`ls SOMATICSNIPER/*$SAMPLE*VEP | grep -v MTOR | grep -v TAU`

STRELKA_SNV=`ls strelka/*$SAMPLE-*/results/pass*snv*VEP | grep -v MTOR | grep -v TAU`
STRELKA_INDEL=`ls strelka/*$SAMPLE-*/results/pass*indel*VEP | grep -v MTOR | grep -v TAU`

#EBCALL=`ls EBCall/EBCall/*$SAMPLE*txt | grep -v MTOR | grep -v TAU`

less $MUTECT | grep -v intron_variant | grep -v ^# | cut -f 1 | sort | uniq | grep -v Uploaded_variation > $NAME.NO_INTRON.EXOMA.MUTECT

less $SOMSNIPER | grep -v intron_variant | grep -v ^# | cut -f 1 | sort | uniq | grep -v Uploaded_variation > $NAME.NO_INTRON.EXOMA.SOMSNIPER

less $STRELKA_SNV | grep -v intron_variant | grep -v ^# | cut -f 1 | sort | uniq | grep -v Uploaded_variation > $SAMPLE.TMP
less $STRELKA_INDEL | grep -v intron_variant | grep -v ^# | cut -f 1 | sort | uniq | grep -v Uploaded_variation > $SAMPLE.TMP2
cat $SAMPLE.TMP $SAMPLE.TMP2 | sort | uniq > $NAME.NO_INTRON.EXOMA.STRELKA
rm $SAMPLE.TMP $SAMPLE.TMP2

#less $EBCALL | grep -v intron_variant | grep -v ^CHR | cut -f 1-4 -d " " | sed 's/chr//g' | awk -F " " '{ print $1 "_" $2 "_" $3 "/" $4}' | sort | uniq >  $NAME.NO_INTRON.EXOMA.EBCALL

#less $EBCALL | grep -v ^CHR | cut -f 1-4 -d " " | sed 's/chr//g' | awk -F " " '{ print $1 "_" $2 "_" $3 "/" $4}' | sort | uniq >  $NAME.EXOMA.EBCALL

done

