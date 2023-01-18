# Convert BAM files to SAM format, which is compatible to HTSeq.
# Inputs: Path to BAM files
# Outputs: SAM files in result folder


# Parse arguments for input path, output path
echo "=====> Converting BAM to SAM"
while getopts 'b:s:' flag
do 
    case "${flag}" in 
        (b) BAM_PATH=${OPTARG};;
        (s) SAM_PATH=${OPTARG};;
        (:) 
            case ${OPTARG} in 
                (b) exit 9999;;
                (s) exit 9999;;
            esac;;
    esac
done


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
    (samtools sort -n -O SAM $f >$SAM_PATH/$ACC_NO.sam) || echo "Err: Cannot convert BAM to SAM $ACC_NO" >>$ENCODE_log/bam2sam.log
}

# Call bam_to_sam() on BAM files folder
for f in $BAM_PATH/*.bam; do
    bam_to_sam $f &
done
wait

echo "=====> Finished mRNA pipeline"
