# All subject tissues ID are stored in $REFGEN/epi_ids_list.main.txt
awk '{ a[NR]=$0 }
       END{ for(i=1;i<=NR;i++)
              for(j=i+1;j<=NR;j++)
		      if (a[i] != a[j]) print a[i], a[j] }' $ENCODE_Escript/epi_ids.txt > all_pairs.txt
mkdir $ENCODE_Escript/pair_chunks
split -l 11 --numeric-suffixes $ENCODE_Escript/all_pairs.txt $ENCODE_Escript/pair_chunks/chunk_
