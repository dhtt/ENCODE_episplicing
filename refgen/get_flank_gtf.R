library("GenomicRanges", quietly = TRUE)
library("rtracklayer", quietly = TRUE)
library("data.table", quietly = TRUE)
library("dplyr", quietly = TRUE)
library("Rsubread", quietly = TRUE)
library("stringr", quietly = TRUE)
library("optparse", quietly = TRUE)

# ==== PREPARATION ====
# Parse arguments for input/output settings
option_list <- list(
  make_option(c("-g", "--reference_genome"),
    type = "character", 
    default = "refgen/reference_genome.gtf",
    help = "path to flattened reference genome", 
    metavar = "character"
  )
  )
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

refgen = import.gff(opt$reference_genome)
refgen = as(refgen, "GRanges")
refgen$exonic_part_number = ifelse(is.na(refgen$exonic_part_number), NA, paste("E", refgen$exonic_part_number, sep=''))

refgen_exon = refgen[refgen$type == "exonic_part"]
refgen.flank.start <- GRanges(
  seqnames = seqnames(refgen_exon),
  ranges = IRanges(start = start(refgen_exon)-200, width = 200*2),
  strand = strand(refgen_exon),
  type = "flank_start",
  gene_id = refgen_exon$gene_id,
  transcript_id = refgen_exon$transcript_id,
  exonic_part_number = refgen_exon$exonic_part_number
)
refgen.flank.end <- GRanges(
  seqnames = seqnames(refgen_exon),
  ranges = IRanges(start = end(refgen_exon)-200, width = 200*2),
  strand = strand(refgen_exon),
  type = "flank_end",
  gene_id = refgen_exon$gene_id,
  transcript_id = refgen_exon$transcript_id,
  exonic_part_number = refgen_exon$exonic_part_number
)
refgen_with_flank = c(refgen.flank.start, refgen.flank.end)
head(refgen_with_flank)
min(start(refgen_with_flank))
start(refgen_with_flank[start(refgen_with_flank) < 0]) = 1
min(start(refgen_with_flank))
refgen_with_flank_sorted = sortSeqlevels(refgen_with_flank)
refgen_with_flank_sorted = sort(refgen_with_flank_sorted)

export(refgen_with_flank_sorted, "refgen/reference_genome.corrected.gtf")
