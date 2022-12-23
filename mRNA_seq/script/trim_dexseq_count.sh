# Filter the original count files to contain only gene that have DEUs occuring
#   in 1-25 pairwise comparisons (out of 171 comparisons). Those are non-ubiquitous
#   DEUs. The genes remained after this filter are subjected to DEU analysis for
#   all of samples using DEXSeq for multiple correction purpose.
# Inputs:
#   DEXSeq_dir: Path to folder containing COUNT_folder. COUNT_folder should
#       contain txt-formatted count files
#   gene_for_correction: txt file containing the genes having DEU in 1-25 pairwise
#       comparisons separated by ','
# Outputs: Count files with filtered genes in txt format in COUNT_folder. Count
#   files are named FILTERED_COUNT_folder/sample1_sample2_count.txt

# Define input and output paths
DEXSeq_dir="$1"
gene_for_correction="$2"
COUNT_folder=$DEXSeq_dir/count
FILTERED_COUNT_folder=$DEXSeq_dir/correction/count
mkdir -p $FILTERED_COUNT_folder

# Read in genes for filtering
corrected_gene=$(cat $gene_for_correction)
corrected_gene=$(sed -e 's:,:\\|:g' <<<$corrected_gene)

# Filter count in original count folder by genes
for x in $COUNT_folder/*.txt; do
    (
        echo $x
        f=${x##*/}
        RES_PATH_temp=$FILTERED_COUNT_folder/temp_$f
        RES_PATH=$FILTERED_COUNT_folder/$f
        grep -w $corrected_gene $x >$RES_PATH_temp
        grep -v '+' $RES_PATH_temp >$RES_PATH
        tail -5 $x >>$RES_PATH
    ) &
done
wait
rm $FILTERED_COUNT_folder/temp_*
echo "===> Finish"
