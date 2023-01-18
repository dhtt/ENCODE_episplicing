1. Define flank regions, discard exon bodies using GRanges 
script: /home/dhthutrang/ENCODE/refgen/get_flank_gtf.R 
input: /home/dhthutrang/ENCODE/refgen/GCF_000001405.39_GRCh38.p13_genomic.gtf 
output: /home/dhthutrang/ENCODE/refgen/reference_genome.fl200.2022.gtf

2. Modify NCBI chromosome nomenclature, filter out unconventional chromosomes 
script: /home/dhthutrang/ENCODE/refgen/prepare_refseq_gtf.sh
input: /home/dhthutrang/ENCODE/refgen/reference_genome.fl200.2022.gtf 
output: /home/dhthutrang/ENCODE/refgen/temp5.gtf

3. Export into gtf format:
script:  /home/dhthutrang/ENCODE/refgen/get_gtf_refgen.R
input: /home/dhthutrang/ENCODE/refgen/temp5.gtf
output: /home/dhthutrang/ENCODE/refgen/hg38.ncbiRefSeq.2022.gtf

3. Flatten gtf:
script: dexseq_prepare_annotation.py 
input: hg38.ncbiRefSeq.2022.gtf 
output: hg38.ncbiRefSeq.2022.flattened.gtf 

