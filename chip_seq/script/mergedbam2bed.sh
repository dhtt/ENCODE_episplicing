# Convert BAM files to BED format
# Inputs: Path to BAM-formatted merged alignment files in alignment_files/merged
# Outputs: Alignments in BED format in alignment_files/bed 


# Define input and output paths
INPUT_PATH=chip_seq/alignment_files/merged
mkdir chip_seq/alignment_files/bed
OUTPUT_PATH=chip_seq/alignment_files/bed

# Convert BAM to BED ===
echo "Convert bam to bed"
for file in $INPUT_PATH/*
do (
    ID=${file%%.*}
    ID=${ID##*/}
    echo $ID
    bedtools bamtobed -i $INPUT_PATH/$ID.bam > $OUTPUT_PATH/$ID.bed
) & 
done
wait
echo "DONE!"
