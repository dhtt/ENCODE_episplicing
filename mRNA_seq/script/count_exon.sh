echo "=====> Begin counting"
cd "$(dirname "$0")"

if (test -f mRNA_seq_BEDaccessionnumber.txt && test -f epi_ids.txt)
then
    echo "...File requirement met."
else
    echo "Err: Check mRNA_seq_BEDaccessionnumber.txt"
    exit 1
fi

# COUNTING =====================================
dexseq_count() {
    f="$p".sam
    ACC_NO=${f%%.*}
    echo "counting $ACC_NO" >> log_dexseqcount.txt
    (python ~/R/x86_64-pc-linux-gnu-library/3.6/DEXSeq/python_scripts/dexseq_count.py $REFGEN/reference_genome.gtf $f backupcount_NCBI/"$ACC_NO"_count.txt) || (echo "Err: Counting $ACC_NO failed" >> log_dexseqcount.txt)
}

while read p; 
do
    dexseq_count &
done < mRNA_seq_BEDaccessionnumber.txt #temp.txt is accession_number (with order sometimes)
wait

echo "=====> Changing names"
# Change name to contain epigenome ID =====================================
for f in backupcount_NCBI/*count.txt
do
    ACC_NO=${f#*/}
    ACC_NO=${ACC_NO%_*}
    ID=$(grep -ih $ACC_NO epi_ids.txt)
    prefix_ID=${ID%%_*}
    NEWNAME="backupcount_NCBI"/"${prefix_ID}"_"$ACC_NO"_"count.txt"
    mv $f $NEWNAME
done

echo "=====> Finished counting"
