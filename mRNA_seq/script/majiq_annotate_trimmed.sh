H=$1
MAJIQ_res='/home/dhthutrang/ENCODE/mRNA_seq/script/majiq_to_be_mapped/*'
REF_GEN=$ENCODE_REFGEN/reference_genome.gtf
MAJIQ_annotated_histone='/home/dhthutrang/ENCODE/mRNA_seq/script/majiq_annotated/'$H
for FILE in $MAJIQ_res
do
    NAME=${FILE##*/}
    NAME=${NAME%%.*}
    OUTPUT=$MAJIQ_annotated_histone'/'$NAME'.bed'
    HIS_ANNOT='/home/dhthutrang/ENCODE/chip_seq/'$H'/flank/'$NAME'.txt.fl.txt'
    echo $H $HIS_ANNOT $OUTPUT 
    bedtools intersect -a $HIS_ANNOT -b $FILE -wo -loj -bed > $OUTPUT
done
wait

echo 'DONE'

