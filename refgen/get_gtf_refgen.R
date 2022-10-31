library(GenomicRanges)
library(rtracklayer)
library(data.table)

setwd("/home/dhthutrang/ENCODE/refgen")
gtf_file = import.gff("temp5.gtf")
head(gtf_file)

#unique(gtf_file$type) # filter only gene, exon/ start codon for promoter

refgen_exon = gtf_file[gtf_file$type %in% c('gene','exon'), c('type', 'gene_id', 'transcript_id', 'exon_number')]
head(refgen_exon)
export(refgen_exon, "hg19.ncbiRefSeq.gtf")
export(refgen_exon, "hg19.ncbiRefSeq.2022.gtf")
