# Perform Differential Exon Usage (DEU) analysis for all pairs of samples using DEXSeq
# Inputs:
#   DEXSEQ_PATH: Path to folder containing COUNT_folder. COUNT_folder should contain txt-formatted count files
#   CHUNK_PATH: Path to folder where the comparisons are splitted into chunks for parallel processing. Chunk files are
#   named chunk_xx.
#   REFGEN_PATH: Path to reference genome .gtf file
# Outputs: DEXSeq results for each pairwise comparison. The results include Rdata/sample1_sample2.RData, 
#   normedcount/sample1_sample2_normedcount.csv, res/sample1_sample2_res.csv. The folders "Rdata", "normedcount",
#    and "res" are generated in the parent directory of COUNT_folder


# Define input path and split workload into chunks for parallel process
echo "===> START DEXSEQ-ING ALL PAIRS"

# Parse arguments for input path, output path
while getopts 'd:c:g:' flag
do 
    case "${flag}" in 
        (d) DEXSEQ_PATH=${OPTARG};;
        (c) CHUNK_PATH=${OPTARG};;
        (g) REFGEN_PATH=${OPTARG};;
        (:) 
            case ${OPTARG} in 
                (d) exit 9999;;
                (c) exit 9999;;
                (g) exit 9999;;
            esac;;
    esac
done

COUNT_folder=$DEXSEQ_PATH/count
RData_folder=$DEXSEQ_PATH/Rdata


#######################################
# Perform DEXSeq for mulitple comparisons in parallel
# Arguments:
#    a tissue pair in each line of the chunk file
# Outputs:
#    DEXSeq results
#######################################
execute_dexseq_parallel() {
    epi1=$(echo $p | awk '{split($0, a, " "); print a[1]}')
    epi2=$(echo $p | awk '{split($0, a, " "); print a[2]}')
    if (ls $COUNT_folder | grep "$epi1") && (ls $COUNT_folder | grep "$epi2"); then
        (
            if (ls $RData_folder | grep "$epi1"_"$epi2"); then
                (
                    echo "Pair "$epi1"_"$epi2" already exists"
                )
            else
                (
                    Rscript mRNA_seq/script/DEXSeq_analysis.R -a "$epi1" -b "$epi2" -f $COUNT_folder -g $REFGEN_PATH -n 8
                    echo "DEXSeq-ing pair "$epi1"_"$epi2""
                )
            fi
        )
    else
        echo "Pair "$epi1"_"$epi2" does not exist"
    fi
}

# Call execute_dexseq_parallel()
for chunk in $CHUNK_PATH/*; do
    while read p; do
        execute_dexseq_parallel &
    done <$chunk
    wait
    echo "===> FINISHED DEXSEQ-ING ALL PAIRS IN "$chunk""
done
echo "End Time: $(date)"
