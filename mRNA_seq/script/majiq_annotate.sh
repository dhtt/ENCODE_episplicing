MAJIQ_res='/home/dhthutrang/ENCODE/mRNA_seq/script/majiq_to_be_mapped/*'
REF_GEN=$ENCODE_REFGEN/reference_genome.gtf
MAJIQ_annotated='/home/dhthutrang/ENCODE/mRNA_seq/script/majiq_annotated/'
for FILE in $MAJIQ_res
do
    NAME=${FILE##*/}
    NAME=${NAME%%.*}
    OUTPUT=$MAJIQ_annotated''$NAME'.bed'
    echo $OUTPUT
    bedtools intersect -a $REF_GEN -b $FILE -wo -loj -bed > $OUTPUT
done
wait

echo 'DONE'
