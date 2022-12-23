# Generate the sample pairs from 19 tissues
# Inputs: none
# Outputs: pair_chunk folder where the comparisons are splitted into chunks
#      for parallel processing. Chunk files are named chunk_xx.
# Globals:
#   ENCODE_EXP: Path to mrna_seq working dir

DEU_script=$ENCODE_EXP/script
epigenome_list=$DEU_script/epi_ids.txt

# All subject tissues ID are stored in $REFGEN/epi_ids.txt
awk '{ a[NR]=$0 }
       END{ for(i=1;i<=NR;i++)
              for(j=i+1;j<=NR;j++)
		      if (a[i] != a[j]) print a[i], a[j] }' $epigenome_list >all_pairs.txt
mkdir $DEU_script/pair_chunks
split -l 11 --numeric-suffixes $DEU_script/all_pairs.txt $DEU_script/pair_chunks/chunk_
