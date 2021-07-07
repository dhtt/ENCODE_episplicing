H=$1
mkdir 'majiq_annotated/'$H
MAJIQ_res='/home/dhthutrang/ENCODE/mRNA_seq/script/majiq_to_be_mapped/*'
REF_GEN=$ENCODE_REFGEN/reference_genome.gtf
MAJIQ_annotated_histone='/home/dhthutrang/ENCODE/mRNA_seq/script/majiq_annotated/'$H
for FILE in $MAJIQ_res
do (
    NAME=${FILE##*/}
    NAME=${NAME%%.*}
    OUTPUT=$MAJIQ_annotated_histone'/'$NAME'.bed'
    OUTPUT_temp=$MAJIQ_annotated_histone'/'$NAME'_temp.bed'
    HIS_ANNOT='/home/dhthutrang/ENCODE/chip_seq/'$H'/flank/bed/'$NAME'.bed'
    echo $H $HIS_ANNOT $OUTPUT 
    bedtools intersect -a $HIS_ANNOT -b $FILE -wo -loj -bed > $OUTPUT_temp
    bedtools groupby -i $OUTPUT_temp -g 1-5 -c 7,8,11 -o collapse,collapse,collapse > $OUTPUT
    )&
done
wait

echo 'DONE'

#cd $MAJIQ_annotated_histone
#mkdir collapsed
#for f in *.bed
#do
#bedtools groupby -i $f -g 1-5 -c 7,8,11 -o collapse,collapse,collapse > collapsed/$f
#done
#echo "Finish"

rm majiq_annotated/$H/*_temp.bed
find ./majiq_annotated/$H -size  0 -print -delete
#wc -l *


