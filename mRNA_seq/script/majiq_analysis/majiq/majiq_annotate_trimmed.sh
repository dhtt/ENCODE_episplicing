#annot_majiq
mkdir majiq_annotated/majiq
all_H=("H3K27ac" "H3K27me3" "H3K36me3" "H3K4me1" "H3K4me3" "H3K9me3")
for H in ${all_H[*]}
do
    mkdir 'majiq_annotated/'$H
done

REF_GEN=$ENCODE_REFGEN/reference_genome.fl200.gtf #to reference_genome.gtf if exon body
MAJIQ_res='/home/dhthutrang/ENCODE/mRNA_seq/script/majiq_to_be_mapped/*'
MAJIQ_annotated='/home/dhthutrang/ENCODE/mRNA_seq/script/majiq_annotated/'
for FILE in $MAJIQ_res
do (
    NAME=${FILE##*/}
    NAME=${NAME%%.*}
    OUTPUT=$MAJIQ_annotated'/majiq/'$NAME'.bed'
    bedtools intersect -a $ENCODE_REFGEN/reference_genome.fl200.gtf -b majiq_to_be_mapped/aorta_H1.deltapsi.tsv -wo -loj -bed | bedtools groupby -g 1-5,7,9 -c 11,12,15 -o collapse,collapse,collapse > $OUTPUT

    for H in ${all_H[*]}
    do (
        HIS_file='/home/dhthutrang/ENCODE/chip_seq/'$H'/flank/bed/'$NAME'.bed'
        awk 'BEGIN { FS=";"} {print $5, $6}' $HIS_file > 'majiq_annotated/'$H'/'$NAME'.bed'
    ) 
    done
)&
done

wait
echo 'DONE'
