# Convert BAM files to SAM format, which is compatible to HTSeq.
# Inputs: Path to BAM files
# Outputs: SAM files in result folder
# Globals:
#     ENCODE_EXP: Path to mrna_seq working dir

# Define input and output paths
echo "=====> Converting BAM to SAM"
BAM_folder=$1
SAM_folder=$ENCODE_EXP/sam_filesÏ

#######################################
# Convert BAM to SAM in parallel
# Arguments:
#    Path to BAM file
# Outputs:
#    SAM files
#######################################
bam_to_sam() {
    ACC_NO=${f%%.*}
    ACC_NO=${ACC_NO##*/}
    echo $ACC_NO
    echo "Converting $ACC_NO" >>$ENCODE_log/bam2sam.log
    (samtools sort -n -O SAM $f >$SAM_folder/$ACC_NO.sam) || echo "Err: Cannot convert BAM to SAM $ACC_NO" >>$ENCODE_log/bam2sam.log
}

# Call bam_to_sam() on BAM files folder
for f in $BAM_folder/*.bam; do
    bam_to_sam $f &
done
wait

echo "=====> Finished mRNA pipeline"
