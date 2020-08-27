cat /dev/null > execute_manorm.log
histone_type=$1
INPUT_PATH=$ENCODE_HIS/$histone_type
echo "Start Time: $(date)" >> execute_manorm.log
echo "===> START MANORM-ING ALL PAIRS IN $histone_type"
mkdir $INPUT_PATH/manorm_result
mkdir $INPUT_PATH/pair_chunks
split -l 11 --numeric-suffixes $INPUT_PATH/all_pairs.txt $INPUT_PATH/pair_chunks/chunk_
ls $INPUT_PATH

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
