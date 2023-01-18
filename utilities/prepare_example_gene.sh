# Gather the raw DEXSeq and MAnorm results for a gene for plotting Figure 4 in utilities/get_example_gene_plot.R
# Inputs:
#   EPI1: first epigenome
#   EPI2: second epigenome
#   GENE_ID: gene which should be plotted
#   REFGEN: Path to the raw reference genome file with all original transcripts (GTF format)
#   REFGEN_DEXSEQ: Path to the DEXSeq-collapsed reference genome file (GTF format)
# Outputs: A folder containing the DEXSeq/MAnorm results for all histone available for this 2 epigenomes


# Parse arguments for input path, output path
echo "===> START DEXSEQ-ING ALL PAIRS"

REFGEN=refgen/hg38.ncbiRefSeq.gtf 
REFGEN_DEXSEQ=refgen/reference_genome.2021.gtf

while getopts 'a:b:g:r:d:' flag
do 
    case "${flag}" in 
        (a) EPI1=${OPTARG};;
        (b) EPI2=${OPTARG};;
        (g) GENE_ID=${OPTARG};;
        (r) REFGEN=${OPTARG};;
        (d) REFGEN_DEXSEQ=${OPTARG};;
        (:) 
            case ${OPTARG} in 
                (a) exit 9999;;
                (b) exit 9999;;
                (g) exit 9999;;
                (r) exit 9999;;
                (d) exit 9999;;
            esac;;
    esac
done

EXAMPLE_GENE_PATH=general_analysis_results/example_gene
mkdir -p $EXAMPLE_GENE_PATH

#Prepare exp file
echo "====> Prepare exp file"
grep "\"$GENE_ID\"" $REFGEN >$EXAMPLE_GENE_PATH/"$GENE_ID"_NCBI.gtf
sed -i 's/_id//g' $EXAMPLE_GENE_PATH/"$GENE_ID"_NCBI.gtf

grep "groupID" mRNA_seq/dexseqcount/res/"$EPI1"_"$EPI2"_res.csv >temp1.txt
grep $GENE_ID$'\t' mRNA_seq/dexseqcount/res/"$EPI1"_"$EPI2"_res.csv >>temp1.txt
grep "\"$GENE_ID\"" $REFGEN_DEXSEQ >temp2.txt
paste -d '\t' temp1.txt temp2.txt >$EXAMPLE_GENE_PATH/"$EPI1"_"$EPI2"_res.csv
rm temp*txt

#Prepare his file
echo "====> Prepare his file"
for his in chip_seq/H*/; do
    (
        temp=${his%*/}
        type=${temp##*/}
        grep "\"$GENE_ID\"" chip_seq/$type/exon_pi/exon_"$EPI1"_"$EPI2".txt >$EXAMPLE_GENE_PATH/exon_"$EPI1"_"$EPI2"."$type".txt
        grep "\"$GENE_ID\"" chip_seq/$type/exon_pi/pi_"$EPI1"_"$EPI2".txt >$EXAMPLE_GENE_PATH/pi_"$EPI1"_"$EPI2"."$type".txt
    )
done

#Finishing
tar -cvf example_"$EPI1"_"$EPI2"_"$GENE_ID".tar.gz $EXAMPLE_GENE_PATH
