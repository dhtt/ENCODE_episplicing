cat /dev/null > execute_dexseq.log
echo "===> START DEXSEQ-ING ALL PAIRS"
f="$1"

execute_dexseq_parallel() {
    epi1=$(echo $p | awk '{split($0, a, " "); print a[1]}')
    epi2=$(echo $p | awk '{split($0, a, " "); print a[2]}')
    if (ls $f | grep "$epi1" ) && (ls $f | grep "$epi2" )
    then (
        if (ls "$f"/RData | grep "$epi1"_"$epi2" )
        then (
            echo "Pair "$epi1"_"$epi2" already exists" >> execute_dexseq.log
        )
        else (
            Rscript $ENCODE_Escript/DEXSeq_analysis.R -a "$epi1" -b "$epi2" -f  "$f" -g $ENCODE_REFGEN/reference_genome.gtf -n 8
            echo "DEXSeq-ing pair "$epi1"_"$epi2"" >> execute_dexseq.log
        )
        fi
        )
    else echo "Pair "$epi1"_"$epi2" does not exist" >> execute_dexseq.log
    fi
}

for chunk in pair_chunks/*
do
    while read p; 
    do execute_dexseq_parallel &
    done < $chunk
    wait
    echo "===> FINISHED DEXSEQ-ING ALL PAIRS IN "$chunk"" >> execute_dexseq.log
done
