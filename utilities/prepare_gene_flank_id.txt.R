library(GenomicRanges)
library(rtracklayer)
library(data.table)
library(dplyr)

#Prepare epi_ids_list.main.txt
file = fread("/Users/dhthutrang/Documents/BIOINFO/Episplicing/ENCODE_episplicing/file_transfer/mrna_submetadata.tsv")
file$tissue = sapply(file$`Biosample term name`, function(x) paste(strsplit(x, split=' ')[[1]], collapse = ''))
temp = file[,c('tissue', 'File accession')]
file_name = as.data.frame(paste(temp$tissue, temp$`File accession`, sep='_'))
fwrite(file_name, "/Users/dhthutrang/Documents/BIOINFO/Episplicing/ENCODE_episplicing/file_transfer/epi_ids.txt",
       col.names = FALSE, row.names = FALSE, quote = FALSE)

unique(temp$tissue)
fwrite(as.data.frame(unique(temp$tissue)), "/Users/dhthutrang/Documents/BIOINFO/Episplicing/ENCODE_episplicing/file_transfer/epi_ids_list.main.txt",
       col.names = FALSE, row.names = FALSE, quote = FALSE)

#Prepare gene_id.txt
refgen.exon = import.gff("/Users/dhthutrang/Documents/BIOINFO/Episplicing/ENCODE_episplicing/refgen/reference_genome.gtf")
refgen.exon = as(refgen.exon, "GRanges")
refgen.exon = refgen.exon[refgen.exon$type == "exonic_part"]
gene_id = unique(refgen.exon$gene_id)
gene_id = as.data.frame(sort(gene_id))
fwrite(gene_id, "/Users/dhthutrang/Documents/BIOINFO/Episplicing/ENCODE_episplicing/utilities/gene_id.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)

#Prepare flank_id.txt
flank_id = as.data.frame(cbind(refgen.exon$gene_id, refgen.exon$exonic_part_number))
flank_id$V2 = paste("E", flank_id$V2, sep='')
fwrite(flank_id, 
       "/Users/dhthutrang/Documents/BIOINFO/Episplicing/ENCODE_episplicing/utilities/flank_id.txt", 
       col.names = FALSE, row.names = FALSE, quote = FALSE, sep='\t')
