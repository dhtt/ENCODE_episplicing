INPUT_PATH=$1
mkdir $ENCODE_HIS/alignment_files/bed
OUTPUT_PATH=$ENCODE_HIS/alignment_files/bed

#=== CONVERT BAM TO BED ===
echo "Convert bam to bed"
for file in $INPUT_PATH/*
do (
    ID=${file%%.*}
    ID=${ID##*/}
    echo $ID
    bedtools bamtobed -i $INPUT_PATH/$ID.bam > $OUTPUT_PATH/$ID.bed
)
done
