# Perform Differential Exon Usage (DEU) analysis for all pairs of samples
#   using DEXSeq
# Inputs:
#   DEXSeq_dir: Path to folder containing COUNT_folder. COUNT_folder should
#       contain txt-formatted count files
#   pair_chunks: Path to folder where the comparisons are splitted into chunks
#       for parallel processing. Chunk files are named chunk_xx.
# Outputs: DEXSeq results for each pairwise comparison. The results include
#   Rdata/sample1_sample2.RData, normedcount/sample1_sample2_normedcount.csv,
#   res/sample1_sample2_res.csv. The result folders "Rdata", "normedcount",
#   and "res" are generated in the parent directory of COUNT_folder
# Globals:
#   ENCODE_EXP: Path to mrna_seq working dir
#   ENCODE_REFGEN: Path to reference genome dir

# Define input path and split workload into chunks for parallel process
cat /dev/null >execute_dexseq.log
echo "===> START DEXSEQ-ING ALL PAIRS"
DEXSeq_dir="$1"
pair_chunks="$2"
REFGEN_path=$ENCODE_REFGEN/reference_genome.gtf
DEXSEQ_script=$ENCODE_EXP/script/DEXSeq_analysis.R
COUNT_folder=$DEXSeq_dir/count
RData_folder=$DEXSeq_dir/Rdata
Rscript $DEXSEQ_script

#######################################
# Perform DEXSeq for mulitple comparisons in parallel
# Globals:
#    pair_chunks
# Arguments:
#    None
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
                    echo "Pair "$epi1"_"$epi2" already exists" >>execute_dexseq.log
                )
            else
                (
                    Rscript $DEXSEQ_script -a "$epi1" -b "$epi2" -f $COUNT_folder -g $REFGEN_path -n 8
                    echo "DEXSeq-ing pair "$epi1"_"$epi2"" >>execute_dexseq.log
                )
            fi
        )
    else
        echo "Pair "$epi1"_"$epi2" does not exist"
    fi
}

# Call execute_dexseq_parallel()
for chunk in $pair_chunks/*; do
    while read p; do
        execute_dexseq_parallel &
    done <$chunk
    wait
    echo "===> FINISHED DEXSEQ-ING ALL PAIRS IN "$chunk"" >>execute_dexseq.log
done
echo "End Time: $(date)" >>$INPUT_PATH/execute_dexseq.log
