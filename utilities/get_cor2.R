library(data.table, quietly=TRUE)
library(dplyr, quietly=TRUE)
library(boot)
library(stats)
library(parallel)
library(tidyverse)
library("doMC")
#setwd("/Users/dhthutrang/Documents/BIOINFO/Episplicing/ENCODE_episplicing/flank")
doMC::registerDoMC(cores = 17)

# histone_type_list = c("H3K27ac", "H3K27me3", "H3K36me3", "H3K4me1", "H3K4me3", "H3K9me3")
histone_type_list = c("H3K36me3")
get_colname <- function(filename_list, option='his'){
  name = sapply(filename_list, function(x) strsplit(x, split='/'))
  name = sapply(name, function(x) x[length(x)][[1]])
  if (option=="his"){
    name = sapply(name, function(x) strsplit(x, split="\\.")[[1]][1])
  }
  else if (option=="exp"){
    name = sapply(name, function(x) strsplit(x, split=".res.csv")[[1]][1])
  }
  name = sapply(name, function(x) paste(sort(strsplit(x, split="_")[[1]] ), collapse = "_"))
  names(name) = NULL
  return(name)
}

# FILTER WITH NEW CORRECTED GENES AND FIRST EXON
filter_genes = function(df, filter_genes_path="combined_df_exon.RDS", filter="deu"){
  # combined_df_exon = readRDS('/Users/trangdo/Documents/Episplicing/ENCODE_episplicing/utilities/combined_df_exon_10perc.RDS')
  combined_df_exon = readRDS(filter_genes_path)
  corrected_genes = unique(unlist(combined_df_exon[combined_df_exon['deu'] == T, "gene_id"]))
  
  filtered_df = as.data.frame(df[df$gene_id %in% corrected_genes, ]) #Filter by corrected genes[only 10 percs]
  combined_df = as.data.frame(combined_df_exon[combined_df_exon$gene_id %in% corrected_genes, ])
  filtered_df = filtered_df[order(filtered_df$gene_id), ] #Order so rownames overlap
  combined_df = combined_df[order(combined_df$gene_id), ]
  filtered_df = filtered_df[combined_df$first_exon == F, ] #Filter out first exons
  combined_df = combined_df[combined_df$first_exon == F, ]
  
  final_filtered_df = filtered_df #Filter by deu/dhm
  if (filter == "deu") final_filtered_df[unlist(combined_df[[filter]] == F), 3:length(final_filtered_df)] = 0
  
  #------Check if filter by DEU work------
  print("========= CHECK FILTER =========")
  print(paste("Gene IDs match:", table(final_filtered_df$gene_id == filtered_df$gene_id)))
  # lapply(list(filtered_df, final_filtered_df), function(x) table(x == 0))
  print(paste("Overlap before after:", table(final_filtered_df == filtered_df)))
  # final_filtered_df[final_filtered_df$adiposetissue_aorta != filtered_df$adiposetissue_aorta, c('gene_id', 'exon_id')]
  # final_filtered_df$adiposetissue_aorta[final_filtered_df$gene_id == 'DLG1' & final_filtered_df$exon_id == 'E015']
  # filtered_df$adiposetissue_aorta[filtered_df$gene_id == 'DLG1' & filtered_df$exon_id == 'E015']
  #------
  
  return(as.data.table(final_filtered_df))
}

#===== PREPARE HIS FILE (6 TOTAL) =====
print("===== PREPARE HIS FILE (6 TOTAL) =====")
his_id = read.csv("flank_id.2021.txt", sep='\t', header = FALSE)
get_all_pairs.his <- function(all_pairs.his, his){
  pair.his_list = vector("list", length(all_pairs.his))
  for (i in 1:length(all_pairs.his)){
    # for (i in 1:1){
    print(paste("Pair: ", i, sep=''))
    pair.his = all_pairs.his[[i]]
    pair.his = fread(pair.his)
    id = as.data.frame(do.call(rbind, lapply(pair.his$V9, function(x) strsplit(x, split='"', fixed=T)[[1]][c(2, 6)])))
    colnames(id) = c('gene', 'exon')
    pair.his = pair.his %>%
      dplyr::mutate(
        gene = id$gene, exon = id$exon, type = V3,
        p_val = as.numeric(as.character(V11)),
        m_val = dplyr::if_else(p_val <= 0.05, 
                               true = abs(as.numeric(as.character(V10))), false = 0)
      ) %>%
      dplyr::select(p_val)
      dplyr::ungroup() 
    if (i == 1) pair.his_id = pair.his[, c('gene', 'exon')]
    print(head(pair.his, 50))
    pair.his_list[[i]] = pair.his
  }
  lapply(pair.his_list, function(x) print(dim(x)))
  pair.his_list = as.data.frame(cbind(pair.his_id, as.data.frame(do.call(cbind, pair.his_list))))
  # saveRDS(pair.his_list, paste('pair.his_list_', his, '.RDS', sep=''))
  return(pair.his_list)
}

get_all_pairs.his_list <- function(histone_type_list){
  all_pairs.his_list = vector("list", length(histone_type_list))
  for (j in 1:length(histone_type_list)){
    his = histone_type_list[[j]]
    all_pairs.his = list.files(paste("/home/dhthutrang/ENCODE/chip_seq", his, "flank/fl", sep='/'), pattern = '.txt', full.names = TRUE)
    print(all_pairs.his)
    colname_his = c("gene_id", "exon_id", get_colname(all_pairs.his, "his")) 
    all_pairs.his.sig = get_all_pairs.his(all_pairs.his, his)
    colnames(all_pairs.his.sig) = colname_his
    
    print(dim(all_pairs.his.sig))
    print(head(all_pairs.his.sig))
    
    all_pairs.his_list[[j]] = as.data.table(all_pairs.his.sig)
  }
  return(all_pairs.his_list)
}
filter_all_his_list <- function(his_list, histone_type_list, filter_genes_path){
  all_filtered_df = vector("list")
  for (i in 1:length(histone_type_list)){
    histone = histone_type_list[i]
    his_df = his_list[[i]]
    print(i)
    print(histone)
    all_filtered_df[[i]] = filter_genes(df = his_df, filter_genes_path = filter_genes_path, filter = histone)
  }
  names(all_filtered_df) = histone_type_list
  return(all_filtered_df)
}

all_pairs.his_list = get_all_pairs.his_list(histone_type_list)
#saveRDS(all_pairs.his_list, "/home/dhthutrang/ENCODE/flank/all_pairs.his_list.RDS")
#all_pairs.his_list_ = readRDS("/home/dhthutrang/ENCODE/flank/all_pairs.his_list.RDS")
# saveRDS(all_pairs.his_list, "/home/dhthutrang/ENCODE/flank/new_df/all_pairs.his_list.RDS")
#all_pairs.his_list_ = lapply(histone_type_list, function(x) readRDS(paste('pair.his_list_', x, '.RDS', sep='')))
#names(all_pairs.his_list_) = histone_type_list

all_pairs.his_list_ = vector("list", length(histone_type_list))
for (j in 1:length(histone_type_list)){
  his = histone_type_list[[j]]
  all_pairs.his = readRDS(paste('pair.his_list_', his, '.RDS', sep=''))
  print(his)
  file.names = list.files(paste("/home/dhthutrang/ENCODE/chip_seq", his, "flank/fl", sep='/'), pattern = '.txt', full.names = TRUE)
  colnames(all_pairs.his) = c("gene_id", "exon_id", get_colname(file.names, "his"))
  all_pairs.his[is.na(all_pairs.his)] = 0
  # if (his == "H3K36me3"){
  #     print('CHECK GSTM3 all_pairs.his$aorta_CD4positivealphabetaTcell')
  #     print(all_pairs.his$aorta_CD4positivealphabetaTcell[all_pairs.his$gene_id =="GSTM3"])
  #     print('CHECK FGFR2 all_pairs.his$aorta_CD4positivealphabetaTcell')
  #     print(all_pairs.his$aorta_CD4positivealphabetaTcell[all_pairs.his$gene_id =="FGFR2"])
  #     }
  all_pairs.his_list_[[j]] = as.data.table(all_pairs.his)
}


#all_pairs.his_list_ = all_pairs.his_list_[c(1,2,3,5,6)] #Leave out H3K4me1
# all_pairs.his_list_flt_10 = filter_all_his_list(his_list = all_pairs.his_list_, histone_type_list=histone_type_list, filter_genes_path="combined_df_exon_10perc.RDS")
# all_pairs.his_list_flt_90 = filter_all_his_list(all_pairs.his_list_, histone_type_list, "combined_df_exon_90perc.RDS")
all_pairs.his_list_flt_90 = filter_all_his_list(all_pairs.his_list_, histone_type_list, "combined_df_exon_90_final.RDS")
# saveRDS(all_pairs.his_list_flt_10, "all_pairs.his_list_flt_10.RDS")
saveRDS(all_pairs.his_list_flt_90, "all_pairs.his_list_flt_90.RDS")
