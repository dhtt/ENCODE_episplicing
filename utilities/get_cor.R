library(data.table, quietly=TRUE)
library(dplyr, quietly=TRUE)
library(boot)
library(stats)
library(parallel)
library("doMC")
#setwd("/Users/dhthutrang/Documents/BIOINFO/Episplicing/ENCODE_episplicing/flank")
doMC::registerDoMC(cores = 17)

histone_type_list = c("H3K27ac", "H3K27me3", "H3K36me3", "H3K4me3", "H3K9me3")
histone_type_list = c("H3K36me3", "H3K27ac")
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

# #===== PREPARE EXP FILE (1 FOR ALL HIS TYPES) =====
# print("===== PREPARE EXP FILE (1 FOR ALL HIS TYPES) =====")
# # all_pairs.exp = list.files("/home/dhthutrang/ENCODE/mRNA_seq/dexseqcount/res", full.names = TRUE)
# get_all_pairs.exp <- function(all_pairs.exp){
#   colname_exp = c("gene_id", "exon_id", get_colname(all_pairs.exp, "exp"))
#   pair.exp_list = vector("list", length(all_pairs.exp))
#   for (i in 1:length(all_pairs.exp)){
#     print(paste("Pair: ",i, sep=''))
#     pair.exp = all_pairs.exp[[i]]
#     #fwrite(exp_id, "/home/dhthutrang/ENCODE/utilities/exp_id.txt", col.names = FALSE, row.names = FALSE, quote = FALSE, sep='\t')
#     pair.exp_list[[i]] = fread(pair.exp)[,c("stat", "padj")]
#   }
#   pair.exp_list = lapply(pair.exp_list,
#                          function(x) {
#                            x = x %>%
#                              mutate(
#                                exp = dplyr::if_else(padj <= 0.1 & !is.na(padj), 
#                                                     true = stat, false = 0.0)) %>%
#                              dplyr::select(exp)
#                          })
#   pair.exp_list = as.data.frame(pair.exp_list)
#   exp_id = fread("/home/dhthutrang/ENCODE/utilities/exp_id.txt", sep = '\t', quote=FALSE, header = FALSE)
#   # print(paste("COMPARE LENGTH", dim(exp_id), dim(pair.exp_list), sep=' '))
#   pair.exp_list = as.data.frame(cbind(exp_id, pair.exp_list))
#   pair.exp_list = pair.exp_list[order(pair.exp_list$V1), ]
#   colnames(pair.exp_list) = colname_exp
#   return(as.data.table(pair.exp_list))
# }
# 
# # all_pairs.exp = get_all_pairs.exp(all_pairs.exp)
# # saveRDS(all_pairs.exp, "/home/dhthutrang/ENCODE/flank/all_pairs.exp.RDS")
# # all_pairs.exp = readRDS("/Users/trangdo/Documents/Episplicing/ENCODE_episplicing/flank/all_pairs.exp.RDS")
# all_pairs.exp_ = readRDS("/home/dhthutrang/ENCODE/flank/all_pairs.exp.RDS")
# # saveRDS(all_pairs.exp, "/home/dhthutrang/ENCODE/flank/new_df/all_pairs.exp.RDS")
# 
# # FILTER WITH NEW CORRECTED GENES AND FIRST EXON
# filter_genes = function(df, filter_genes_path="combined_df_exon.RDS", filter="deu"){
#   # combined_df_exon = readRDS('/Users/trangdo/Documents/Episplicing/ENCODE_episplicing/utilities/combined_df_exon_10perc.RDS')
#   combined_df_exon = readRDS(filter_genes_path)
#   corrected_genes = unique(unlist(combined_df_exon[combined_df_exon['deu'] == T, "gene_id"]))
#   
#   filtered_df = as.data.frame(df[df$gene_id %in% corrected_genes, ]) #Filter by corrected genes[only 10 percs]
#   combined_df = as.data.frame(combined_df_exon[combined_df_exon$gene_id %in% corrected_genes, ])
#   filtered_df = filtered_df[order(filtered_df$gene_id), ] #Order so rownames overlap
#   combined_df = combined_df[order(combined_df$gene_id), ]
#   filtered_df = filtered_df[combined_df$first_exon == F, ] #Filter out first exons
#   combined_df = combined_df[combined_df$first_exon == F, ]
#   
#   final_filtered_df = filtered_df #Filter by deu/dhm
#   final_filtered_df[unlist(combined_df[[filter]] == F), 3:length(final_filtered_df)] = 0
#   
#   #------Check if filter by DEU work------
#   print("========= CHECK FILTER =========")
#   print(paste("Gene IDs match:", table(final_filtered_df$gene_id == filtered_df$gene_id)))
#   # lapply(list(filtered_df, final_filtered_df), function(x) table(x == 0))
#   print(paste("Overlap before after:", table(final_filtered_df == filtered_df)))
#   # final_filtered_df[final_filtered_df$adiposetissue_aorta != filtered_df$adiposetissue_aorta, c('gene_id', 'exon_id')]
#   # final_filtered_df$adiposetissue_aorta[final_filtered_df$gene_id == 'DLG1' & final_filtered_df$exon_id == 'E015']
#   # filtered_df$adiposetissue_aorta[filtered_df$gene_id == 'DLG1' & filtered_df$exon_id == 'E015']
#   #------
#   
#   return(as.data.table(final_filtered_df))
# }
# 
# # all_pairs.exp_flt_10 = filter_genes(all_pairs.exp_, filter_genes_path="combined_df_exon_10perc.RDS", filter="deu")
# # all_pairs.exp_flt_90 = filter_genes(all_pairs.exp_, filter_genes_path="combined_df_exon_90perc.RDS", filter="deu")
# all_pairs.exp_flt_90 = filter_genes(all_pairs.exp_, filter_genes_path="combined_df_exon_90_final.RDS", filter="deu")
# # saveRDS(all_pairs.exp_flt_10, "all_pairs.exp_flt_10.RDS")
# # saveRDS(all_pairs.exp_flt_90, "all_pairs.exp_flt_90.RDS")
# all_pairs.exp_flt_10 = readRDS("all_pairs.exp_flt_10.RDS")
# # all_pairs.exp_flt_90 = readRDS("all_pairs.exp_flt_90.RDS")

#===== PREPARE HIS FILE (6 TOTAL) =====
print("===== PREPARE HIS FILE (6 TOTAL) =====")
his_id = read.csv("flank_id.2021.txt", sep='\t', header = FALSE)
get_all_pairs.his <- function(all_pairs.his){
  pair.his_list = vector("list", length(all_pairs.his))
  for (i in 1:length(all_pairs.his)){
    # for (i in 1:1){
    print(paste("Pair: ", i, sep=''))
    pair.his = all_pairs.his[[i]]
    pair.his = fread(pair.his)
    if (i==1) {
      his_id = as.data.frame(lapply(pair.his$V9, function(x) strsplit(x, split='"', fixed=T)[[1]][c(2, 6)]))
      colnames(his_id) = c("V1", "V2")
      print(head(his_id))
      }
    pair.his = pair.his %>%
      mutate(
        p_val = as.numeric(as.character(V11)),
        m_val = dplyr::if_else(p_val <= 0.05, 
                               true = abs(as.numeric(as.character(V10))), false = 0)
      ) %>%
      dplyr::select(m_val)
    pair.his_list[[i]] = pair.his
  }
  pair.his_list = as.data.table(pair.his_list)
  pair.his_list = pair.his_list %>%
    group_by(group = gl(n()/2, 2)) %>%
    summarise_all(max) %>%
    dplyr::select(-group)
  pair.his_list = as.data.frame(cbind(his_id, pair.his_list))
  pair.his_list = pair.his_list[order(pair.his_list$V1),]
  return(pair.his_list)
}

get_all_pairs.his_list <- function(histone_type_list){
  all_pairs.his_list = vector("list", length(histone_type_list))
  for (j in 1:length(histone_type_list)){
    # for (j in 1:1){
    his = histone_type_list[[j]]
    all_pairs.his = list.files(paste("/home/dhthutrang/ENCODE/chip_seq", his, "flank/fl", sep='/'), pattern = '.txt', full.names = TRUE)
    print(all_pairs.his)
    colname_his = c("gene_id", "exon_id", get_colname(all_pairs.his, "his")) 
    all_pairs.his.sig = get_all_pairs.his(all_pairs.his)
    colnames(all_pairs.his.sig) = colname_his
    print(all_pairs.his.sig$aorta_CD8positivealphabetaTcell[all_pairs.his.sig$gene_id == "FGFR2"])
    
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
saveRDS(all_pairs.his_list, "/home/dhthutrang/ENCODE/flank/all_pairs.his_list.RDS")
all_pairs.his_list_ = readRDS("/home/dhthutrang/ENCODE/flank/all_pairs.his_list.RDS")
# saveRDS(all_pairs.his_list, "/home/dhthutrang/ENCODE/flank/new_df/all_pairs.his_list.RDS")

all_pairs.his_list_ = all_pairs.his_list_[c(1,2,3,5,6)] #Leave out H3K4me1
# all_pairs.his_list_flt_10 = filter_all_his_list(his_list = all_pairs.his_list_, histone_type_list=histone_type_list, filter_genes_path="combined_df_exon_10perc.RDS")
# all_pairs.his_list_flt_90 = filter_all_his_list(all_pairs.his_list_, histone_type_list, "combined_df_exon_90perc.RDS")
all_pairs.his_list_flt_90 = filter_all_his_list(all_pairs.his_list_, histone_type_list, "combined_df_exon_90_final.RDS")
# saveRDS(all_pairs.his_list_flt_10, "all_pairs.his_list_flt_10.RDS")
saveRDS(all_pairs.his_list_flt_90, "all_pairs.his_list_flt_90.RDS")
# all_pairs.his_list_flt_10 = readRDS("all_pairs.his_list_flt_10.RDS")
# all_pairs.his_list_flt_90 = readRDS("all_pairs.his_list_flt_90.RDS")

# all_pairs.his_list = readRDS("all_pairs.his_list.RDS")
# head(all_pairs.his_list[[1]], 50)
# 
# his_ids = as.data.frame(read.csv("flank_id.txt", header= FALSE, sep='\t'))
# his_ids = his_ids[order(his_ids$V1, his_ids$V2),]
# all_pairs.exp = as.data.frame(all_pairs.exp)
# all_pairs.exp = all_pairs.exp[order(all_pairs.exp$gene_id, all_pairs.exp$exon_id),]
# table(all_pairs.exp$gene_id ==his_ids$V1)
# 
# #===== CORRELATION WITH RANDOMIZATION =====
# # ------------ Execute analysis ------------
# # all_genes = fread("gene_id.txt", header = FALSE)
# # all_genes = unique(all_pairs.exp_flt$gene_id)
# p_value_calculator <- function(r, nrow){
#   P <- r*sqrt(nrow-2)/sqrt(1-r*r)
#   P <- 2*pt(-abs(P), nrow-2)
#   return(P)
# }
# pearcor_p <- function(exp, his){
#   if (length(unique(exp)) > 1 & length(unique(his)) > 1){
#     p_val = p_value_calculator(cor(exp, his, method = "pearson"), nrow = length(exp))
#     # p_val = p.adjust(p_val, method = "fdr", n=ncol(all_genes))
#     return(p_val)
#   }
#   else {
#     return(NA)
#   }
# }
# 
# 
# pearcor_r <- function(exp, his, n_points){
#   df = as.data.frame(cbind(exp, his))
#   n_sep_point = nrow(unique(df))
#   if (0 %in% apply(df, 1, unique)) n_sep_point = n_sep_point - 1
#   if (n_sep_point >= 3 & length(unique(exp)) > 1 & length(unique(his)) >= n_points){
#     r_val = cor(exp, his, method = "pearson")
#     return(r_val)
#   }
#   else {
#     return(NA)
#   }
# }
# 
# analyze_array <- function(all_pairs.exp, all_pairs.his, option = "p", n_points, all_genes){
#   all_res_pair = vector("list", ncol(all_pairs.exp) - 2)
#   subset_name = colnames(all_pairs.his)
#   print(subset_name)
#   colnames(all_pairs.exp) = gsub('trophoblastcell', 'trophoblast', colnames(all_pairs.exp))
#   all_pairs.exp_subset = all_pairs.exp[, ..subset_name]
#   
#   #for (i in 1:n_pairs){
#   all_res_pair <- foreach( i=1:(length(subset_name)-2), .combine='c', .packages=c('dplyr') ) %dopar% {
#     print(paste("Pair: ", i, sep=''))
#     exp = all_pairs.exp_subset[[i+2]]
#     his = all_pairs.his[[i+2]]
#     data_table = as.data.table(cbind(exp, his))
#     head(data_table)
#     
#     if (option == "p"){
#       res_table = data_table %>%
#         dplyr::group_by(all_pairs.exp$gene_id) %>%
#         dplyr::summarise(res = pearcor_p(exp, his)) %>%
#         dplyr::select(res)
#     }
#     else if (option == "r") {
#       res_table = data_table %>%
#         dplyr::group_by(all_pairs.exp$gene_id) %>%
#         dplyr::summarise(res = pearcor_r(exp, his, n_points)) %>%
#         dplyr::select(res)
#     }
#     #all_res_pair[[i]] = res_table
#   }
#   all_res_pair = as.data.table(all_res_pair)
#   all_res_pair = cbind(all_genes, all_res_pair)
#   colnames(all_res_pair) = c('gene_id', subset_name[3:length(subset_name)])
#   print(head(all_res_pair))
#   return(as.data.frame(all_res_pair))
# }
# analyze_array_list <- function(all_pairs.exp, all_pairs.his_list, method = "p", n_points=2){
#   all_res_list = vector("list", length(histone_type_list)-1 )
#   all_genes = unique(all_pairs.exp$gene_id)
#   for (j in 1:length(histone_type_list)){
#     print(paste("Histone: ", histone_type_list[[j]], sep = ''))
#     all_pairs.his = all_pairs.his_list[[j]]
#     if (method == "p") {
#       all_res_pair = analyze_array(all_pairs.exp, all_pairs.his, option = "p", all_genes=all_genes)
#     }
#     else if (method == "r") {
#       all_res_pair = analyze_array(all_pairs.exp, all_pairs.his, option = "r", n_points=n_points, all_genes=all_genes)
#     }
#     all_res_list[[j]] = all_res_pair
#   }
#   return(all_res_list)
# }
# 
# print("Pearsons-p correlation")
# # all_res_list.pearcor_p = analyze_array_list(all_pairs.exp_flt_10, all_pairs.his_list_flt_10, method = "p")
# # saveRDS(all_res_list.pearcor_p, "/home/dhthutrang/ENCODE/flank/new_df/all_res_list.pearcor_p_10perc.RDS")
# # 
# # all_res_list.pearcor_r = analyze_array_list(all_pairs.exp_flt_10, all_pairs.his_list_flt_10, method = "r")
# # saveRDS(all_res_list.pearcor_r, "/home/dhthutrang/ENCODE/flank/new_df/all_res_list.pearcor_r_10perc.RDS")
# 
# all_res_list.pearcor_p = analyze_array_list(all_pairs.exp_flt_90, all_pairs.his_list_flt_90, method = "p")
# saveRDS(all_res_list.pearcor_p, "/home/dhthutrang/ENCODE/flank/new_df/all_res_list.pearcor_p_90_final.RDS")
# 
# all_res_list.pearcor_r = analyze_array_list(all_pairs.exp_flt_90, all_pairs.his_list_flt_90, method = "r")
# saveRDS(all_res_list.pearcor_r, "/home/dhthutrang/ENCODE/flank/new_df/all_res_list.pearcor_r_90_final.RDS")
# # 
# # all_res_list.pearcor_p = analyze_array_list(all_pairs.exp, all_pairs.his_list, method="r", n_points=5)
# # saveRDS(all_res_list.pearcor_p, "/home/dhthutrang/ENCODE/flank/new_df/all_res_list.pearcor_r5.RDS")
# # 
# # all_res_list.pearcor_p = analyze_array_list(all_pairs.exp, all_pairs.his_list, method="r", n_points=10)
# # saveRDS(all_res_list.pearcor_p, "/home/dhthutrang/ENCODE/flank/new_df/all_res_list.pearcor_r10.RDS")
# # 
# # all_res_list.pearcor_p = analyze_array_list(all_pairs.exp, all_pairs.his_list, method="r", n_points=15)
# # saveRDS(all_res_list.pearcor_p, "/home/dhthutrang/ENCODE/flank/new_df/all_res_list.pearcor_r15.RDS")
# 
# 
# # temp = fread("/Users/dhthutrang/Documents/BIOINFO/Episplicing/ENCODE_episplicing/utilities/CD4positivealphabetaTcell_endodermalcell.txt.fl.txt")
# # temp[grep("STIM1", temp$V9), `V10`]
# # name = temp$V9[seq(1, nrow(temp), 2)]
# # 
# # ids = lapply(name, function(x){
# #   split_name = str_split(x, '\"')[[1]]
# #   gene_id = split_name[2]
# #   exon_id = split_name[length(split_name) -1]
# #   return(c(gene_id, exon_id))
# # })
# # tail(ids)
# # his_id2 = as.data.frame(do.call("rbind", ids))
# # head(his_id2)
# # 
# # temp1 = temp %>%
# #   dplyr::select(V10) %>%
# #   group_by(group = gl(n()/2, 2)) %>%
# #   summarise_all(max) %>%
# #   dplyr::select(-group)
# # head(temp3)
# # his_id2 = read.csv("flank_id.txt", sep='\t', header = FALSE)
# # temp2 = cbind(temp1, his_id)
# # temp3 = cbind(temp1, his_id2)
# # temp2[temp2$V1 == "STIM1", 1]
# # temp3[temp3$V1 == "STIM1", 1]
# # 
# # fwrite(his_id2, "flank_id.txt", quote = FALSE, col.names = FALSE, sep='\t')
# 
# 
