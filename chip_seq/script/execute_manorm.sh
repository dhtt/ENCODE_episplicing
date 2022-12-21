# Execute MAnorm to compare the cells/tissues in a pairwise manner for each histone type
# Inputs: The files used for the pairwise comparison are defined in histone_type/pair_chunks/. 
#    These include the peak files in peak_files/merged/ and alignment files in 
#    alignment_files/bed, both in BED format.
# Outputs: For each pairwise comparison, a folder containing MAnorm results are stored in 
#    histone_type/manorm_result/pairwise_comparison/
# Globals:
#     ENCODE_HIS: Path to chip_seq working dir


# Define input path and split workload into chunks for parallel process
cat /dev/null > execute_manorm.log
histone_type=$1
mkdir $ENCODE_HIS/$histone_type
INPUT_PATH=$ENCODE_HIS/$histone_type

echo "Start Time: $(date)" >> execute_manorm.log
echo "===> START MANORM-ING ALL PAIRS IN $histone_type"
mkdir $INPUT_PATH/manorm_result
mkdir $INPUT_PATH/pair_chunks
split -l 11 --numeric-suffixes $INPUT_PATH/all_pairs.txt $INPUT_PATH/pair_chunks/chunk_


#######################################
# Perform MAnorm for mulitple comparisons in parallel
# Globals:
#    INPUT_PATH
# Arguments:
#    None
# Outputs:
#    MAnorm results
#######################################
execute_manorm_parallel() {
    peak_file1=$(echo $p | awk '{split($0, a, " "); print a[1]}')
    read_file1=$(echo $p | awk '{split($0, a, " "); print a[2]}')
    peak_file2=$(echo $p | awk '{split($0, a, " "); print a[3]}')
    read_file2=$(echo $p | awk '{split($0, a, " "); print a[4]}')
    epi1=$(echo $p | awk '{split($0, a, " "); print a[5]}')
    epi2=$(echo $p | awk '{split($0, a, " "); print a[6]}')
    if (ls $INPUT_PATH/manorm_result | grep "$epi1"_"$epi2" )
    then (
        echo "Pair "$epi1"_"$epi2" already exists" >> execute_manorm.log
        )
    else (
        echo "MAnorm-ing pair "$epi1"_"$epi2": "$peak_file1"|"$peak_file2"|"$read_file1"|"$read_file2"" >> $INPUT_PATH/execute_manorm.log
        manorm --p1 "$peak_file1" --p2 "$peak_file2" --r1 "$read_file1" --r2 "$read_file2" -o $INPUT_PATH/manorm_result/"$epi1"_"$epi2"
        )
    fi
}


# Call execute_manorm_parallel()
for chunk in $INPUT_PATH/pair_chunks/*
do
    while read p;
    do
    execute_manorm_parallel &
    done < $chunk
    wait
    echo "===> FINISHED manorm-ING ALL PAIRS IN "$chunk"" >> $INPUT_PATH/execute_manorm.log
done
echo "End Time: $(date)" >> $INPUT_PATH/execute_manorm.log
