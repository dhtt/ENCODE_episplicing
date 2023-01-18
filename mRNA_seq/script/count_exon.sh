# Generate count files from SAM-formatted read files using DEXSeq/HTSeq script.
# Inputs: Path to SAM_folder of SAM-formatted read files
# Outputs: Count files in txt format in COUNT_folder. Count files are named in 
#   COUNT_folder/sample1_sample2_count.txt
# Globals:
#   ENCODE_EXP: Path to mrna_seq working dir

# Parse arguments for input path, output path
while getopts 'i:o:g:' flag
do 
    case "${flag}" in 
        (i) INPUT_PATH=${OPTARG};;
        (o) OUTPUT_PATH=${OPTARG};;
        (g) REFGEN_PATH=${OPTARG};;
        (:) 
            case ${OPTARG} in 
                (i) exit 9999;;
                (o) exit 9999;;
                (g) exit 9999;;
            esac;;
    esac
done

# Start exon reads counting
echo "generate_exon_count.sh: STARTED"

dexseq_count() {
    ID=${file%%.*}
    ID=${ID##*/}
    OUTPUT_FILENAME=$OUTPUT_PATH/"$ID"_count.txt
    echo "Counting $ID" >> DEU_log.txt
    if [ -e "$OUTPUT_FILENAME" ]
    then 
        echo "Already counted $ID" >> DEU_log.txt    
    else
        (python ~/R/x86_64-pc-linux-gnu-library/3.6/DEXSeq/python_scripts/dexseq_count.py $REFGEN_PATH $file $OUTPUT_FILENAME) || (echo "Err: Counting $ID failed" >> DEU_log.txt)
    fi
}
for file in $INPUT_PATH/*.sam
do
    dexseq_count & 
done
wait
echo "generate_exon_count.sh: DONE"