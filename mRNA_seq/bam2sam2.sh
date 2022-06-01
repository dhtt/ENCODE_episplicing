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

#for f in $BAM_folder/*.bam
#do
#    bam_to_sam $f 
#done
#wait

samtools sort -n -O SAM /home/dhthutrang/ENCODE/mRNA_seq/raw_data/6/ENCFF081DYZ.bam > $SAM_folder/ENCFF081DYZ.sam
samtools sort -n -O SAM /home/dhthutrang/ENCODE/mRNA_seq/raw_data/6/ENCFF268XDH.bam > $SAM_folder/ENCFF268XDH.sam
samtools sort -n -O SAM /home/dhthutrang/ENCODE/mRNA_seq/raw_data/6/ENCFF359YOQ.bam > $SAM_folder/ENCFF359YOQ.sam
samtools sort -n -O SAM /home/dhthutrang/ENCODE/mRNA_seq/raw_data/6/ENCFF588TQX.bam > $SAM_folder/ENCFF588TQX.bam
#bam_to_sam /home/dhthutrang/ENCODE/mRNA_seq/raw_data/6/ENCFF085TWC.bam
#bam_to_sam /home/dhthutrang/ENCODE/mRNA_seq/raw_data/6/ENCFF717MVQ.bam

echo "=====> Finished mRNA pipeline"
