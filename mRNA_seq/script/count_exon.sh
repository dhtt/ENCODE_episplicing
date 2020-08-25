echo "=====> Begin counting"
cd "$(dirname "$0")"
mkdir $ENCODE_EXP/dexseqcount #dexseqcount instead of backupcount_NCBI

# COUNTING =====================================
dexseq_count() {
    ID=${file%%.*}
    ID=${ID##*/}
    echo $ID
    echo $file
    echo "counting $ID" >> log_dexseqcount.txt
    # (python ~/R/x86_64-pc-linux-gnu-library/3.6/DEXSeq/python_scripts/dexseq_count.py $ENCODE_REFGEN/reference_genome.gtf $file $ENCODE_EXP/dexseqcount/"$ID"_count.txt) || (echo "Err: Counting $ID failed" >> log_dexseqcount.txt)
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
    # ID=$(grep -ih $ACC_NO epi_ids.txt)
    # prefix_ID=${ID%%_*}
    # NEWNAME="backupcount_NCBI"/"${prefix_ID}"_"$ACC_NO"_"count.txt"
    # mv $f $NEWNAME
done
echo "=====> Finished counting"
