library(GenomicRanges)
library(rtracklayer)
library(data.table)
library(dplyr)

refgen.exon = import.gff("/Users/dhthutrang/Documents/BIOINFO/Episplicing/files/gene_level/reference_genome.exon.gtf")
refgen.exon = as(refgen.exon, "GRanges")
start(refgen.exon)
strand(refgen.exon)
refgen.exon$transcripts

refgen.flank.start <- GRanges(
  seqnames = seqnames(refgen.exon),
  ranges = IRanges(start = start(refgen.exon)-100, width = 100*2),
  strand = strand(refgen.exon),
  type = "flank_start",
  gene = refgen.exon$gene_id,
  transcript = refgen.exon$transcripts,
  exon = paste("E", refgen.exon$exonic_part_number, sep='')
)
refgen.flank.end <- GRanges(
  seqnames = seqnames(refgen.exon),
  ranges = IRanges(start = end(refgen.exon)-100, width = 100*2),
  strand = strand(refgen.exon),
  type = "flank_end",
  gene = refgen.exon$gene_id,
  transcript = refgen.exon$transcripts,
  exon = paste("E", refgen.exon$exonic_part_number, sep='')
)
refgen.flank = c(refgen.flank.start, refgen.flank.end)
min(start(refgen.flank))
start(refgen.flank[start(refgen.flank) < 0]) = 1
temp = sortSeqlevels(refgen.flank)
temp = sort(temp)
# temp = as.data.frame(temp)

export(temp, "reference_genome.fl100.gtf")








