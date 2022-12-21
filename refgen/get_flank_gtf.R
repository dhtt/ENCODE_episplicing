library(GenomicRanges)
library(rtracklayer)
library(data.table)
library(dplyr)
library(Rsubread)
library(stringr)


# refgen = import.gff("/home/dhthutrang/ENCODE/refgen/reference_genome.2021_.gtf")

# refgen = as(refgen, "GRanges")
# refgen$exonic_part_number = ifelse(is.na(refgen$exonic_part_number), NA, paste("E", refgen$exonic_part_number, sep=''))
# # start(refgen.exon)
# # strand(refgen.exon)
# # refgen.exon$transcript_id
# refgen_exon = refgen[refgen$type == "exonic_part"]
# refgen.flank.start <- GRanges(
#   seqnames = seqnames(refgen_exon),
#   ranges = IRanges(start = start(refgen_exon)-200, width = 200*2),
#   strand = strand(refgen_exon),
#   type = "flank_start",
#   gene_id = refgen_exon$gene_id,
#   transcript_id = refgen_exon$transcript_id,
#   exonic_part_number = refgen_exon$exonic_part_number
# )
# refgen.flank.end <- GRanges(
#   seqnames = seqnames(refgen_exon),
#   ranges = IRanges(start = end(refgen_exon)-200, width = 200*2),
#   strand = strand(refgen_exon),
#   type = "flank_end",
#   gene_id = refgen_exon$gene_id,
#   transcript_id = refgen_exon$transcript_id,
#   exonic_part_number = refgen_exon$exonic_part_number
# )
# refgen_with_flank = c(refgen.flank.start, refgen.flank.end)
# head(refgen_with_flank)
# min(start(refgen_with_flank))
# start(refgen_with_flank[start(refgen_with_flank) < 0]) = 1
# min(start(refgen_with_flank))
# temp = sortSeqlevels(refgen_with_flank)
# temp = sort(temp)

# export(temp, "/home/dhthutrang/ENCODE/refgen/reference_genome.2021.corrected.gtf")


# refgen = import.gff("/home/dhthutrang/ENCODE/refgen/reference_genome.2021.corrected.gtf")
# paste(unique(refgen$gene_id), collapse = "\\|")

# refgen = import.gff("/home/dhthutrang/ENCODE/refgen/GCF_000001405.39_GRCh38.p13_genomic.gtf")
# refgen_exon = refgen[refgen$type == "exon" & refgen$gbkey == "mRNA"]
# new_flank_start <- GRanges(
#   seqnames = seqnames(refgen_exon),
#   ranges = IRanges(start = start(refgen_exon)-200, width = 200*2),
#   strand = strand(refgen_exon),
#   type = "exon",
#   source = refgen_exon$source,
#   gene_id = refgen_exon$gene_id,
#   transcript_id = refgen_exon$transcript_id, 
#   exon_number = refgen_exon$exon_number
# )
# new_flank_end <- GRanges(
#   seqnames = seqnames(refgen_exon),
#   ranges = IRanges(start = end(refgen_exon)-200, width = 200*2),
#   strand = strand(refgen_exon),
#   type = "exon",
#   source = refgen_exon$source,
#   gene_id = refgen_exon$gene_id,
#   transcript_id = refgen_exon$transcript_id, 
#   exon_number = refgen_exon$exon_number
# )
# refgen_with_flank = c(new_flank_start, new_flank_end)
# refgen_with_flank = sortSeqlevels(refgen_with_flank)
# refgen_with_flank = sort(refgen_with_flank)
# start(refgen_with_flank[start(refgen_with_flank) < 0]) = 1

# export(refgen_with_flank, "/home/dhthutrang/ENCODE/refgen/reference_genome.fl200.2022.gtf")


# refgen = import.gff("/home/dhthutrang/ENCODE/refgen/hg19.ncbiRefSeq.2022.flattened.gtf")
# refgen_exon = refgen[refgen$type == "exonic_part"]
# refgen_exon_alphabet = refgen_exon[order(refgen_exon$gene_id, refgen_exon$exonic_part_number)] %>% as.data.frame()
# result_list = vector("list")
# result_paths =  list.files("/home/dhthutrang/Krebs/Reddy/data/mRNA_seq/DEXSeq_output/csv", full.names = TRUE, pattern = 'res.csv')
# for (i in seq(length(result_paths))){
#   result = fread(file)
#   result = result[order(result$groupID, result$featureID),]
#   result$stat = ifelse(result$padj <= 0.05, result$stat, 0) 
#   result_list[[i]] = result$stat

# head(refgen_exon_alphabet)
# refgen_exon_alphabet = do.call('cbind', list(refgen_exon_alphabet, as.data.frame(result_list)))
# colnames(refgen_exon_alphabet)[13 : ncol(refgen_exon_alphabet)] = sapply(result_paths, function(x){
#   y = strsplit(x, '/', fixed = TRUE)[[1]]
#   return(y[length(y)])
#   })
# refgen_exon = as(refgen_exon_alphabet, "GRanges")
# refgen_exon = sortSeqlevels(refgen_exon)
# refgen_exon = sort(refgen_exon)
# refgen_exon = as.data.frame(refgen_exon)
# head(refgen_exon[, c(1, 2, 3, 13:ncol(refgen_exon))])
# fwrite(refgen_exon[, c(1, 2, 3, 13:ncol(refgen_exon))], col.names = FALSE, row.names = FALSE, sep = '\t', na = 0, quote = FALSE,
#   "/home/dhthutrang/Krebs/Reddy/data/mRNA_seq/dexseq_output_annotated/dexseq_res.bed")
# # bedtools intersect -a /home/dhthutrang/ENCODE/refgen/reference_genome.fl200.2021_.corrected.gtf -b /home/dhthutrang/Krebs/Reddy/data/mRNA_seq/dexseq_output_annotated/dexseq_res.bed -wo -loj -bed > /home/dhthutrang/Krebs/Reddy/data/mRNA_seq/dexseq_output_annotated/dexseq_output_annotated.bed
# #  bedtools groupby -i /home/dhthutrang/Krebs/Reddy/data/mRNA_seq/dexseq_output_annotated/dexseq_output_annotated.bed -g 1-9 -c 13 -o max > /home/dhthutrang/Krebs/Reddy/data/mRNA_seq/dexseq_output_annotated/dexseq_output_annotated.grouped.bed 


refgen <- import.gff('/home/dhthutrang/Krebs/episplice-pipeline/utils/prepare_refgen/hg19.ncbiRefSeq.2022.gtf')
id_lookup_table <- unique(as.data.frame(refgen)[, c('gene_id', 'transcript_id')])

flatten_refgen <- flattenGTF('/home/dhthutrang/Krebs/episplice-pipeline/utils/prepare_refgen/hg19.ncbiRefSeq.2022.gtf', 
  GTF.featureType='exon', GTF.attrType = "transcript_id", method="chop")
saveRDS(flatten_refgen, "flatten_refgen.RDS")

flatten_refgen_collapse <- flatten_refgen  %>%
  dplyr::group_by(Chr, Start, End, Strand) %>%
  dplyr::mutate(
    gene_id = id_lookup_table$gene_id[match(GeneID, id_lookup_table$transcript_id)],
    transcript_id = paste(GeneID, collapse = '+'),
    type = 'exonic_part') %>%
  dplyr::select(-GeneID) %>%
  unique()  %>%
  ungroup() %>%
  dplyr::group_by(gene_id, transcript_id)  %>%
  dplyr::mutate(
    exon_number = 1:length(transcript_id)
    ) %>%
  ungroup() 
saveRDS(flatten_refgen_collapse, "flatten_refgen_collapse.RDS")
flatten_refgen_collapse = readRDS("flatten_refgen_collapse.RDS")
flatten_refgen_collapse$type = 'exonic_part'

aggregate_genes <- flatten_refgen_collapse %>%
  dplyr::group_by(gene_id) %>%
  dplyr::summarise(
    Chr = Chr, Start = min(Start), End = max(End),
    Strand = Strand, transcript_id=NA, exon_number=NA, 
    type="aggregate_gene"
    ) %>%
  ungroup() %>%
  unique() %>%
  dplyr::select(Chr, Start, End, Strand, gene_id, transcript_id, type, exon_number)
saveRDS(aggregate_genes, "aggregate_genes.RDS")
aggregate_genes = readRDS("aggregate_genes.RDS")

refgen_new = as.data.frame(rbind(aggregate_genes, flatten_refgen_collapse))
colnames(refgen_new) = c('chr', 'start', 'end', 'strand', 'gene_id', 'transcripts', 'type', 'exonic_part_number')
refgen_new$exonic_part_number = str_pad(refgen_new$exonic_part_number, 3, pad="0")
refgen_new = sortSeqlevels(as(refgen_new, "GRanges"))
refgen_new = sort(refgen_new)
export(refgen_new, "/home/dhthutrang/Krebs/episplice-pipeline/utils/prepare_refgen/refgen_collapsed.gtf")



# test <- flattenGTF("/home/dhthutrang/Krebs/episplice-pipeline/utils/prepare_refgen/test.collapsed_transcript.gtf",
#   GTF.featureType='exon', GTF.attrType = "gene_id", method="chop"
# )
# test_refgen = import.gff("/home/dhthutrang/Krebs/episplice-pipeline/utils/prepare_refgen/test.collapsed_transcript.gtf")
# id_lookup_table <- unique(as.data.frame(test_refgen)[, c('gene_id', 'transcript_id')])

# test_collapse <- test  %>%
#   dplyr::group_by(Chr, Start, End, Strand) %>%
#   dplyr::mutate(
#     gene_id = paste(GeneID, collapse = '+'),
#     transcript_id = id_lookup_table$transcript_id[match(GeneID, id_lookup_table$gene_id)],
#     type = 'exon') %>%
#   dplyr::select(-GeneID) %>%
#   unique()  %>%
#   ungroup() %>%
#   dplyr::group_by(gene_id, transcript_id)  %>%
#   dplyr::mutate(
#     exon_number = 1:length(transcript_id)
#     )
# print(as.data.frame(test_collapse))
# as(flatten_refgen_collapse, "GRanges")
