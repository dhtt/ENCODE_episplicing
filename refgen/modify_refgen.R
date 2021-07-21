########################################################
## Remove dup genes, overlapping genes, report number ##
########################################################
library(GenomicRanges)
library(rtracklayer)
library(data.table, quietly=TRUE)
library(tidyverse)
library(dplyr, quietly=TRUE)
library(plyr)
library("ggVennDiagram")

refgen.exon = import.gff("/Users/dhthutrang/Documents/BIOINFO/Episplicing/ENCODE_episplicing/refgen/reference_genome.gtf")
refgen.exon = as(refgen.exon, "GRanges")

#==== Gene statistics ====
head(refgen.exon)
# Initial no of genes
length(unique(refgen.exon$gene_id)) #19240
# Overlapping genes
overlapping_genes = unique(refgen.exon$gene_id[grep('+', refgen.exon$gene_id, fixed = TRUE)]) #275
# Dupplicated genes
dup_genes = unique(refgen.exon$gene_id[grep('_', refgen.exon$gene_id, fixed = TRUE)]) #17
new_refgen = refgen.exon[refgen.exon$gene_id %in% c(overlapping_genes, dup_genes) == F,]
# Single exon genes
new_refgen_df = data.frame(new_refgen)
new_refgen_df_multiexon = new_refgen_df %>%
  dplyr::group_by(gene_id) %>%
  dplyr::mutate(n_exon = dplyr::n()-1) %>%
  dplyr::filter(n_exon > 1)

length(unique(new_refgen_df$gene_id)) - length(unique(new_refgen_df_multiexon$gene_id)) #1050 
setdiff(unique(new_refgen_df$gene_id), unique(new_refgen_df_multiexon$gene_id))
#Before: 19240 - After: 17900
  

#==== Exon/Transcript statistics ====
# Exon data new
exon_data = new_refgen_df_multiexon %>%
  dplyr::select(gene_id, n_exon) %>%
  dplyr::mutate(type = 'Processed') %>%
  unique()

# Exon data raw
refgen.raw = import.gff("/Users/dhthutrang/Documents/BIOINFO/Episplicing/ENCODE_episplicing/refgen/GCF_000001405.39_GRCh38.p13_genomic.gtf")
refgen.raw = as(refgen.raw, "GRanges")
refgen.raw = refgen.raw[refgen.raw$type %in% c('exon'),]
refgen.raw_df = data.frame(refgen.raw)
exon_data.raw = refgen.raw_df %>%
  dplyr::select(gene_id, transcript_id) %>%
  dplyr::group_by(gene_id) %>%
  dplyr::mutate(n_exon = dplyr::n(),
                n_transcript = length(unique(transcript_id))) %>%
  dplyr::select(gene_id, n_exon, n_transcript) %>%
  unique() %>%
  dplyr::mutate(n_exon = round(n_exon/n_transcript, 0),
                type = 'Raw') %>%
  dplyr::select(gene_id, n_exon, type) 



#First exon
temp = new_refgen_df_multiexon
all_first_exon = c()
for (gene in unique(temp$gene_id)){
  print(gene)
  temp1 = temp[temp$gene_id == gene, ]
  transcript_ = c()
  transcript_list = lapply(temp1$transcript_id, function(x) str_split(x, '\\+')[[1]])
  strand_ = as.character(unique(temp1$strand))
  n_ = length(temp1$transcript_id)
  first_exon = vector('list', n_)
  
  if (strand_ == '+')
    order_ = seq(1, n_, 1)
  else
    order_ = rev(seq(1, n_, 1))
  
  for (i in order_){
    if (all(transcript_list[[i]] %in% transcript_))
      first_exon[i] = FALSE
    else {
      transcript_ = union(transcript_, unlist(transcript_list[[i]], recursive = T))
      first_exon[i] = TRUE
    }
  }
  
  all_first_exon = c(all_first_exon, first_exon)
}
temp$first_exon = all_first_exon
table(unlist(all_first_exon)) #False: 233694  True: 76099 -> Rate: 24,56%
# export(temp, "reference_genome.2021.gtf")
exon_data.no_firstexon = temp %>%
  dplyr::filter(first_exon == FALSE) %>%
  dplyr::select(gene_id, n_exon) %>%
  dplyr::group_by(gene_id) %>%
  dplyr::mutate(n_exon = n(),
                type='w/o First exon') %>%
  unique()

all_exon_data = rbind(exon_data, exon_data.raw, exon_data.no_firstexon)

plot_1 = ggplot(data=all_exon_data, aes(x = n_exon, fill = type)) +
  geom_histogram(aes(y=(..density..)), binwidth = 1, position="dodge", boundary=0, alpha = 0.75) +
  geom_density(alpha=0.45)+
  xlab("Number of exon bins per gene") + ylab("Occurences") +
  theme_bw() + 
  theme(aspect.ratio=1, plot.margin	= unit(c(0.2,0,3.15,0), "cm"),
        legend.title = element_blank(), legend.position = c(0.8, 0.8),
        legend.box.background = element_rect(color = 'lightgrey')) +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 10)) +
  scale_fill_manual(values=c("#18bc81", "#bddf26", "#6b1bc3")) 
plot_1





