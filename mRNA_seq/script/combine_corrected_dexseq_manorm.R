## ---------------------------
##
## Script name: DEXSeq_analysis.R
##
## Purpose of script: Perform pairwise Differential Exon Usage analysis using DEXSeq R package
##
## Author: Trang Do
##
## Copyright (c) Trang Do, 2018
## Email: dhttrang@bioinformatik.uni-saarland.de
##
## ---------------------------
##
## Notes: This script was used for pairwise DEU analysis. For multiple correction, 
## DEXSeq_analysis_combined.R was used
##
## ---------------------------


library(data.table)
library(dplyr)

setwd("/Users/trangdo/Documents/Episplicing/ENCODE_episplicing/mRNA_seq/script/")
histone_types = c("H3K27ac", "H3K27me3", "H3K36me3", "H3K4me3", "H3K9me3")
#======COMBINE H3K27ac with dexseq ====
#-----LOAD FILES-----
# Histone
his_bed_files = list.files(path="/Users/trangdo/Documents/Episplicing/ENCODE_episplicing/chip_seq/manorm2/sig_peaks/annotated", full.names = T)

combine_annotated_his_peaks = function(his_list){
  his_df_list_ = lapply(his_list, fread)
  his_array_ = lapply(his_df_list_[2:length(his_list)], function(x) return(x[['V10']]))
  his_df = cbind(his_df_list_[[1]], do.call(cbind, his_array_))
  colnames(his_df) = c('chr', 'source', 'type', 'start', 'end', 'score', 'strand', 'quality', 'metadata', histone_types)
  his_df = his_df[his_df$type != "aggregate_gene",]
  
  metadata_df = lapply(his_df$metadata, function(x) unlist(strsplit(x, '"'))[c(2, 6, 8)])
  metadata_df = as.data.frame(do.call("rbind", metadata_df))
  colnames(metadata_df) = c('gene_id', 'exon_id', 'first_exon')
  metadata_df$exon_id[grep('E', metadata_df$exon_id, invert = T)] = paste('E', metadata_df$exon_id[grep('E', metadata_df$exon_id, invert = T)], sep='')
  metadata_df$id = paste(metadata_df$gene_id, metadata_df$exon_id, sep=':')
  his_df = cbind(his_df, metadata_df)
  his_df = his_df %>%
    dplyr::select(-metadata) %>%
    dplyr::group_by(id) %>%
    mutate(first_exon = first_exon[!is.na(first_exon)]) %>%
    dplyr::ungroup() %>%
    mutate(H3K27ac = dplyr::if_else(H3K27ac == '.', FALSE, TRUE),
           H3K27me3 = dplyr::if_else(H3K27me3 == '.', FALSE, TRUE),
           H3K36me3 = dplyr::if_else(H3K36me3 == '.', FALSE, TRUE),
           H3K4me3 = dplyr::if_else(H3K4me3 == '.', FALSE, TRUE),
           H3K9me3 = dplyr::if_else(H3K9me3 == '.', FALSE, TRUE))
  his_df$dhm = Reduce('|', his_df[histone_types])
  return(his_df)
}
his_df = combine_annotated_his_peaks(his_bed_files)
# Exon
corrected_deu = readRDS("multiple_correction_sig_exon_75low.RDS")

#Combined info
combined_df = his_df
combined_df$deu = sapply(unlist(his_df$id), function(x) x %in% corrected_deu)
# saveRDS(combined_df, "combined_df.RDS") #still have first exon 

combined_df_exon = combined_df[combined_df$type == "exonic_part",]
# saveRDS(combined_df_exon, "combined_df_exon.RDS") #still have first exon but only exon

his_df_nofirst = combined_df[combined_df$first_exon == FALSE, ]
combined_df_final = his_df_nofirst
combined_df_final[c(histone_types, 'dhm', 'deu')] = lapply(c(histone_types, 'dhm', 'deu'), function(x) as.factor(combined_df_final[x][[1]]))
# saveRDS(combined_df_final, "combined_df_final.RDS") #no first exon 

#======PEFORM TESTING======
#Fisher test for each histone, all histones, ALL genes 
fisher_all = lapply(c(histone_types, 'dhm'), function(x) 
  return(fisher.test(as.factor(combined_df_final[x][[1]]), as.factor(combined_df_final['deu'][[1]])))
  )
names(fisher_all) = c(histone_types, 'dhm')

#Fisher test for each histone, all histones, EACH genes 
combined_df_final_flank = combined_df_final[combined_df_final$type %in% c('flank_start', 'flank_end'), ]
combined_df_final_exon = combined_df_final[combined_df_final$type == 'exonic_part', ]

get_fisher_res = function(meta_df){
  res = meta_df %>%
    group_by(gene_id) %>%
    summarise(H3K27ac = dplyr::if_else(length(unique(H3K27ac)) == 2 && length(unique(deu)) == 2,
                                       true = fisher.test(H3K27ac, deu)$p.value,
                                       false = 1.0),
              H3K27me3 = dplyr::if_else(length(unique(H3K27me3)) == 2 && length(unique(deu)) == 2,
                                        true = fisher.test(H3K27me3, deu)$p.value,
                                        false = 1.0),
              H3K36me3 = dplyr::if_else(length(unique(H3K36me3)) == 2 && length(unique(deu)) == 2,
                                        true = fisher.test(H3K36me3, deu)$p.value,
                                        false = 1.0),
              H3K4me3 = dplyr::if_else(length(unique(H3K4me3)) == 2 && length(unique(deu)) == 2,
                                       true = fisher.test(H3K4me3, deu)$p.value,
                                       false = 1.0),
              H3K9me3 = dplyr::if_else(length(unique(H3K9me3)) == 2 && length(unique(deu)) == 2,
                                       true = fisher.test(H3K9me3, deu)$p.value,
                                       false = 1.0),
              dhm = dplyr::if_else(length(unique(dhm)) == 2 && length(unique(deu)) == 2,
                                   true = fisher.test(dhm, deu)$p.value,
                                   false = 1.0)) %>%
    as.data.frame()
  rownames(res) = res$gene_id 
  res = select(res, -gene_id)
  
  print('Number of sig genes (Fisher exact test):')
  apply(res, 2, function(x) print(table(x < 0.05)))
  
  res.adj = as.data.frame(apply(res, 2, function(x) p.adjust(x, method = 'fdr')), row.names = rownames(res))
  print('Number of sig genes (Fisher exact test adjusted):')
  apply(res.adj, 2, function(x) print(table(x < 0.05)))
  
  print('All adjusted sig genes (Fisher exact test):')
  paste(unlist(apply(res.adj, 2, function(x) rownames(res.adj)[x <= 0.05])), collapse = ',')
  return(res.adj)
}
fisher_res = get_fisher_res(combined_df_final)
fisher_res_flank = get_fisher_res(combined_df_final_flank)
fisher_res_exon = get_fisher_res(combined_df_final_exon)


