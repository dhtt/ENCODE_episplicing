grep 'gbkey "mRNA"' GCF_000001405.39_GRCh38.p13_genomic.gtf > temp.gtf
awk -F'\t' -vOFS='\t' '{ gsub("NC_0000", "chr", $1); print }' temp.gtf > temp1.gtf
awk -F'\t' -vOFS='\t' '{ gsub("chr0", "chr", $1); print }' temp1.gtf > temp2.gtf
awk -F'\t' -vOFS='\t' '{ gsub("\\.10|\\.11|\\.12|\\.9|\\.14|", "", $1); print }' temp2.gtf > temp3.gtf
awk -F'\t' -vOFS='\t' '{ gsub("chr23", "chrX", $1); gsub("chr24", "chrY", $1); print}' temp3.gtf > temp4.gtf
grep "^chr" temp4.gtf > temp5.gtf
Rscript get_gtf_refgen.R
python dexseq_prepare_annotation.py hg19.ncbiRefSeq.gtf hg19.ncbiRefSeq.DEXSeq.gtf #or flattened instead of DEXSeq
