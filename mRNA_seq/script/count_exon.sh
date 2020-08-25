echo "=====> Begin counting"
cd "$(dirname "$0")"
mkdir $ENCODE_EXP/dexseqcount #dexseqcount instead of backupcount_NCBI

# COUNTING =====================================
dexseq_count() {
    ID=${file%%.*}
    ID=${ID##*/}
    echo $ID
    echo $file
    echo "Counting $ID" >> log_dexseqcount.txt
    if [ -e "$ENCODE_EXP/dexseqcount/"$ID"_count.txt" ]
    then 
        echo "ALready counted $ID" >> log_dexseqcount.txt    
    else
        (python ~/R/x86_64-pc-linux-gnu-library/3.6/DEXSeq/python_scripts/dexseq_count.py $ENCODE_REFGEN/reference_genome.gtf $file $ENCODE_EXP/dexseqcount/"$ID"_count.txt) || (echo "Err: Counting $ID failed" >> log_dexseqcount.txt)
    fi
}
for file in ls $ENCODE_EXP/sam_files/*.sam
do
    dexseq_count &
done
wait 

echo "=====> Changing names"


# Change name to contain epigenome ID =====================================
for f in $ENCODE_EXP/dexseqcount/*count.txt
do
    ACC_NO=${f#*/}
    ACC_NO=${ACC_NO%_*}
    ID=$(grep -ih $ACC_NO $ENCODE_EXP/epi_ids.txt)
    echo $ID
    # prefix_ID=${ID%%_*}
    # NEWNAME="backupcount_NCBI"/"${prefix_ID}"_"$ACC_NO"_"count.txt"
    # mv $f $NEWNAME
done
echo "=====> Finished counting"
