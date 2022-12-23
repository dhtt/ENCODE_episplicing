# Generate count files from SAM-formatted read files using DEXSeq/HTSeq script.
# Inputs: Path to SAM_folder of SAM-formatted read files
# Outputs: Count files in txt format in COUNT_folder. Count files are named in 
#   COUNT_folder/sample1_sample2_count.txt
# Globals:
#   ENCODE_EXP: Path to mrna_seq working dir

# Define input and output paths
echo "=====> Begin counting"
cd "$(dirname "$0")"
SAM_folder=$1
COUNT_folder=$ENCODE_EXP/dexseqcount

#######################################
# Generate count files from SAM files in parallel
# Arguments:
#    Path to SAM file
# Outputs:
#    Count files
#######################################
dexseq_count() {
    ID=${file%%.*}
    ID=${ID##*/}
    file_path=$COUNT_folder/$ID"_count.txt"
    echo "Counting $ID" >>log_dexseqcount.txt
    if [ -e $file_path ]; then
        echo "Already counted $ID" >>log_dexseqcount.txt
    else
        echo $file_path
        (python ~/R/x86_64-pc-linux-gnu-library/3.6/DEXSeq/python_scripts/dexseq_count.py $ENCODE_REFGEN/reference_genome.gtf $file $ENCODE_EXP/dexseqcount/"$ID"_count.txt) || (echo "Err: Counting $ID failed" >>log_dexseqcount.txt)
    fi
}

# Call dexseq_count() on SAM files folder
for file in ls $SAM_folder/*.sam; do
    dexseq_count &
done
wait
echo "=====> Finished counting"
