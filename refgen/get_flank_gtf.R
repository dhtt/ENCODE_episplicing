library(GenomicRanges)
library(rtracklayer)
library(data.table)
library(dplyr)

refgen.exon = import.gff("/Users/dhthutrang/Documents/BIOINFO/Episplicing/ENCODE_episplicing/refgen/reference_genome.gtf")
refgen.exon = as(refgen.exon, "GRanges")
start(refgen.exon)
strand(refgen.exon)
refgen.exon$transcript_id
refgen.exon = refgen.exon[refgen.exon$type == "exonic_part"]

refgen.flank.start <- GRanges(
  seqnames = seqnames(refgen.exon),
  ranges = IRanges(start = start(refgen.exon)-200, width = 200*2),
  strand = strand(refgen.exon),
  type = "flank_start",
  gene_id = refgen.exon$gene_id,
  transcript_id = refgen.exon$transcript_id,
  exonic_part_number = paste("E", refgen.exon$exonic_part_number, sep='')
)
refgen.flank.end <- GRanges(
  seqnames = seqnames(refgen.exon),
  ranges = IRanges(start = end(refgen.exon)-200, width = 200*2),
  strand = strand(refgen.exon),
  type = "flank_end",
  gene_id = refgen.exon$gene_id,
  transcript_id = refgen.exon$transcript_id,
  exonic_part_number = paste("E", refgen.exon$exonic_part_number, sep='')
)
refgen.flank = c(refgen.flank.start, refgen.flank.end)
head(refgen.flank)
min(start(refgen.flank))
start(refgen.flank[start(refgen.flank) < 0]) = 1
min(start(refgen.flank))
temp = sortSeqlevels(refgen.flank)
temp = sort(temp)

export(temp, "reference_genome.fl200.gtf")

file = fread("/Users/dhthutrang/Documents/BIOINFO/Episplicing/ENCODE_episplicing/file_transfer/mrna_submetadata.tsv")
file$tissue = sapply(file$`Biosample term name`, function(x) paste(strsplit(x, split=' ')[[1]], collapse = ''))
temp = file[,c('tissue', 'File accession')]
file_name = as.data.frame(paste(temp$tissue, temp$`File accession`, sep='_'))
fwrite(file_name, "/Users/dhthutrang/Documents/BIOINFO/Episplicing/ENCODE_episplicing/file_transfer/epi_ids.txt",
       col.names = FALSE, row.names = FALSE, quote = FALSE)

unique(temp$tissue)
fwrite(as.data.frame(unique(temp$tissue)), "/Users/dhthutrang/Documents/BIOINFO/Episplicing/ENCODE_episplicing/file_transfer/epi_ids_list.main.txt",
       col.names = FALSE, row.names = FALSE, quote = FALSE)







