echo "=====> Converting BAM to SAM"
BAM_folder=$1
SAM_folder=$ENCODE_EXP/sam_files
echo $BAM_folder
echo $SAM_folder


bam_to_sam(){
    ACC_NO=${f%%.*}
    ACC_NO=${ACC_NO##*/}
    echo $ACC_NO
    echo "Converting $ACC_NO" >> $ENCODE_log/bam2sam.log
    (samtools sort -n -O SAM $f > $SAM_folder/$ACC_NO.sam) || echo "Err: Cannot convert BAM to SAM $ACC_NO" >> $ENCODE_log/bam2sam.log
}

for f in $BAM_folder/*.bam
do
    bam_to_sam $f &
done
wait

echo "=====> Finished mRNA pipeline"
