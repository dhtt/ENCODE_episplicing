# Annotate M-Values from MAnorm results to flank regions in reference genome
# Inputs: xls-formatted MAnorm results stored in histone_type/manorm_result/pairwise_comparison/
# Outputs: fl.txt in histone_type/flank/
# Globals:
#     ENCODE_HIS: Path to chip_seq working dir


# Write log file
cat /dev/null > annotate_manorm.log
echo "Start Time: $(date)" >> annotate_manorm.log
histonetype=$1
HISTONE_PATH=chip_seq/$histonetype


#Get chr, start, end, M-value, p-value, normedcount1, normedcount2
echo "===> Begin annotating"
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


# Annotate all xls.bed in normalizedcounts into $FLANK/histonetype
echo "===> 2: Annotate all xls.bed in normalizedcounts into FLANK/histonetype"
mkdir $HISTONE_PATH/flank
mkdir $HISTONE_PATH/flank
retrieve_epiid(){
    FILE1=${FILE##*/}
    epi1=${FILE1%%_*}
    FILE2=${FILE1#_*}
    epi2=${FILE2#*_}
    epi2=${epi2%%_*}
}
annotate_manorm_parallel_flank() {
    retrieve_epiid
    echo $FILE
    echo $FILE2
    echo $epi1
    echo $epi2
    if (ls $HISTONE_PATH/normalizedcounts| grep "$epi1" ) && (ls $HISTONE_PATH/normalizedcounts | grep "$epi2" )
    then (
        if (ls $HISTONE_PATH/flank | grep "$epi1"_"$epi2")
        then (
            echo "Pair "$epi1"_"$epi2" already annotated" >> annotate_manorm.log
        )
        else (
            echo "Annotating pair $epi1 and $epi2" >> annotate_manorm.log

            REF_GEN_EXON=$ENCODE_REFGEN/reference_genome.fl200.gtf
            HIS_COUNT=$FILE
            ANNOT_HIS_COUNT=$HISTONE_PATH/flank/"$epi1"_"$epi2".txt
            bedtools intersect -a $REF_GEN_EXON -b $HIS_COUNT -wo -loj -bed > $ANNOT_HIS_COUNT
            )
        fi
        )
    else echo "Pair "$epi1"_"$epi2" does not exist" >> annotate_manorm.log
    fi
}

for FILE in $HISTONE_PATH/normalizedcounts/*
do
    annotate_manorm_parallel_flank &
done
echo "End Time: $(date)" >> annotate_manorm.log
wait


# Collapse counts
echo "===> 3: Collapsing annotated counts"
for f in $HISTONE_PATH/flank/*
do (
    echo $f 
    bedtools groupby -i $f -g 1-9 -c 13,14 -o max > $f.fl.txt
) &
done
wait
echo "Finish"
