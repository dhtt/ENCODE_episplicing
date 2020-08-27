cat /dev/null > annotate_manorm.log
echo "Start Time: $(date)" >> annotate_manorm.log
histonetype=$1
HISTONE_PATH=$ENCODE_HIS/$histonetype

echo "===> Begin annotating"
#Get chr, start, end, M-value, p-value, normedcount1, normedcount2
echo "===> 1: Get chr, start, end, M-value, p-value, normedcount1, normedcount2"
mkdir $HISTONE_PATH/normalizedcounts
for x in $HISTONE_PATH/manorm_result/*/*.xls
do (
    echo $x
    f=${x##*/}
    echo $f
    awk '{print $1"\t"$2"\t"$3"\t"$5"\t"$7"\t"$9"\t"$10}' $x > $HISTONE_PATH/normalizedcounts/$f.bed
) &
done
wait
echo "===> Finish"

# #annotate all xls.bed in normalizedcounts into $FLANK/histonetype
# echo "===> 2: Annotate all xls.bed in normalizedcounts into FLANK/histonetype"
# mkdir $HISTONE_PATH/flank
# annotate_manorm_parallel_flank() {
#     epi1=$(echo $p | awk '{split($0, a, " "); print a[5]}')
#     epi2=$(echo $p | awk '{split($0, a, " "); print a[6]}')
#     if (ls $HISTONE_PATH/normalizedcounts| grep "$epi1" ) && (ls $HISTONE_PATH/normalizedcounts | grep "$epi2" )
#     then (
#         if (ls $HISTONE_PATH/flank | grep "exon_"$epi1"_"$epi2"")
#         then (
#             echo "Pair "$epi1"_"$epi2" already annotated" >> annotate_manorm.log
#         )
#         else (
#             echo "Annotating pair "$epi1"_"$epi2" - EXON" >> annotate_manorm.log
#             REF_GEN_EXON=$ENCODE_REFGEN/reference_genome.fl200.gtf
#             HIS_COUNT=$HISTONE_PATH/normalizedcounts/"$epi1"_"$epi2"_all_MAvalues.xls.bed
#             ANNOT_HIS_COUNT=$HISTONE_PATH/flank/"$epi1"_"$epi2".txt
#             bedtools intersect -a $REF_GEN_EXON -b $HIS_COUNT -wo -loj -bed > $ANNOT_HIS_COUNT
#             )
#         fi
#         )
#     else echo "Pair "$epi1"_"$epi2" does not exist" >> annotate_manorm.log
#     fi
# }

# for chunk in $HISTONE_PATH/pair_chunks/*
# do
#     while read p;
#     do
#     annotate_manorm_parallel_flank &
#     done < $chunk
#     wait
#     echo "===> FINISHED ANNOTATING ALL PAIRS IN "$chunk" - EXON" >> annotate_manorm.log
# done
# echo "End Time: $(date)" >> annotate_manorm.log

# for f in $FLANK/$histonetype/*
# do (
#     echo $f 
#     bedtools groupby -i $f -g 1-9 -c 13,14 -o max > $f.fl.txt
# ) &
# done
# wait
# echo "Finish"

