cat /dev/null > execute_manorm.log
echo "Start Time: $(date)" >> execute_manorm.log
echo "===> START MANORM-ING ALL PAIRS"
mkdir manorm_result
mkdir pair_chunks
split -l 11 --numeric-suffixes all_pairs.txt pair_chunks/chunk_

execute_manorm_parallel() {
    peak_file1=$(echo $p | awk '{split($0, a, " "); print a[1]}')
    read_file1=$(echo $p | awk '{split($0, a, " "); print a[2]}')
    peak_file2=$(echo $p | awk '{split($0, a, " "); print a[3]}')
    read_file2=$(echo $p | awk '{split($0, a, " "); print a[4]}')
    epi1=$(echo $p | awk '{split($0, a, " "); print a[5]}')
    epi2=$(echo $p | awk '{split($0, a, " "); print a[6]}')
    if (ls manorm_result | grep "$epi1"_"$epi2" )
    then (
        echo "Pair "$epi1"_"$epi2" already exists" >> execute_manorm.log
        )
    else (
        echo "MAnorm-ing pair "$epi1"_"$epi2": "$peak_file1"|"$peak_file2"|"$read_file1"|"$read_file2"" >> execute_manorm.log
        manorm --p1 "$peak_file1" --p2 "$peak_file2" --r1 "$read_file1" --r2 "$read_file2" -o manorm_result/"$epi1"_"$epi2"
        )
    fi
}

for chunk in pair_chunks/*
do
    while read p;
    do
    execute_manorm_parallel &
    done < $chunk
    wait
    echo "===> FINISHED manorm-ING ALL PAIRS IN "$chunk"" >> execute_manorm.log
done
echo "End Time: $(date)" >> execute_manorm.log

