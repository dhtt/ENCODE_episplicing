library(data.table, quietly=TRUE)
library(dplyr, quietly=TRUE)
library(boot)
library(stats)
library(parallel)
library("doMC")
#setwd("/Users/dhthutrang/Documents/BIOINFO/Episplicing/ENCODE_episplicing/flank")
doMC::registerDoMC(cores = 17)
histone_type_list = c("H3K27ac", "H3K27me3", "H3K36me3", "H3K4me3", "H3K9me3")
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

#===== PREPARE EXP FILE (1 FOR ALL HIS TYPES) =====
print("===== PREPARE EXP FILE (1 FOR ALL HIS TYPES) =====")
all_pairs.exp_ = readRDS("/home/dhthutrang/ENCODE/flank/all_pairs.exp.RDS")

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
  final_filtered_df[unlist(combined_df[[filter]] == F), 3:length(final_filtered_df)] = 0
  
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

all_pairs.exp_flt_90 = filter_genes(all_pairs.exp_, filter_genes_path="combined_df_exon_90_final.RDS", filter="deu")
get_odd_df = function(df){
  df_odd = as.data.frame(apply(df, 2, function(x) x = ifelse(abs(as.numeric(x)) > 0, T, F)))
  df_odd$gene_id = df$gene_id
  df_odd$exon_id = df$exon_id
  return(df_odd)
}
all_pairs.exp_flt_90_odd = get_odd_df(all_pairs.exp_flt_90)
saveRDS(all_pairs.exp_flt_90_odd, "all_pairs.exp_flt_90_odd.RDS")
all_pairs.exp_flt_90_odd = readRDS("all_pairs.exp_flt_90_odd.RDS")

#===== PREPARE HIS FILE (6 TOTAL) =====
print("===== PREPARE HIS FILE (6 TOTAL) =====")
his_id = read.csv("flank_id.txt", sep='\t', header = FALSE)
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
all_pairs.his_list_ = readRDS("/home/dhthutrang/ENCODE/flank/all_pairs.his_list.RDS")

all_pairs.his_list_ = all_pairs.his_list_[c(1,2,3,5,6)] #Leave out H3K4me1
all_pairs.his_list_flt_90 = filter_all_his_list(all_pairs.his_list_, histone_type_list, "combined_df_exon_90_final.RDS")
all_pairs.his_list_flt_90_odd = lapply(all_pairs.his_list_flt_90, function(x) get_odd_df(x))
saveRDS(all_pairs.his_list_flt_90_odd, "all_pairs.his_list_flt_90_odd.RDS")
all_pairs.his_list_flt_90_odd = readRDS("all_pairs.his_list_flt_90_odd.RDS")

#===== CORRELATION WITH RANDOMIZATION =====
# ------------ Execute analysis ------------
p_value_calculator <- function(r, nrow){
  P <- r*sqrt(nrow-2)/sqrt(1-r*r)
  P <- 2*pt(-abs(P), nrow-2)
  return(P)
}
pearcor_p <- function(exp, his){
  if (length(unique(exp)) > 1 & length(unique(his)) > 1){
    p_val = questionr::odds.ratio(table(exp, his))$p
    return(p_val)
  }
  else {
    return(NA)
  }
}
pearcor_r <- function(exp, his, n_points){
  df = as.data.frame(cbind(exp, his))
  if (length(unique(exp)) > 1 & length(unique(his)) > 1){
    r_val = questionr::odds.ratio(table(exp, his))$OR
    return(r_val)
  }
  else {
    return(NA)
  }
}

analyze_array <- function(all_pairs.exp, all_pairs.his, option = "p", n_points, all_genes){
  all_res_pair = vector("list", ncol(all_pairs.exp) - 2)
  subset_name = colnames(all_pairs.his)
  print(subset_name)
  colnames(all_pairs.exp) = gsub('trophoblastcell', 'trophoblast', colnames(all_pairs.exp))
  print(head(all_pairs.exp))
  print(head(all_pairs.his))
  print(head(all_pairs.exp[, subset_name]))
  all_pairs.exp_subset = all_pairs.exp[, subset_name]
  
  #for (i in 1:n_pairs){
  all_res_pair <- foreach( i=1:(length(subset_name)-2), .combine='c', .packages=c('dplyr') ) %dopar% {
    print(paste("Pair: ", i, sep=''))
    exp = all_pairs.exp_subset[[i+2]]
    his = all_pairs.his[[i+2]]
    data_table = as.data.table(cbind(exp, his))
    head(data_table)
    
    if (option == "p"){
      res_table = data_table %>%
        dplyr::group_by(all_pairs.exp$gene_id) %>%
        dplyr::summarise(res = pearcor_p(exp, his)) %>%
        dplyr::select(res)
    }
    else if (option == "r") {
      res_table = data_table %>%
        dplyr::group_by(all_pairs.exp$gene_id) %>%
        dplyr::summarise(res = pearcor_r(exp, his, n_points)) %>%
        dplyr::select(res)
    }
    #all_res_pair[[i]] = res_table
  }
  all_res_pair = as.data.table(all_res_pair)
  all_res_pair = cbind(all_genes, all_res_pair)
  colnames(all_res_pair) = c('gene_id', subset_name[3:length(subset_name)])
  print(head(all_res_pair))
  return(as.data.frame(all_res_pair))
}
analyze_array_list <- function(all_pairs.exp, all_pairs.his_list, method = "p", n_points=2){
  all_res_list = vector("list", length(histone_type_list)-1 )
  all_genes = unique(all_pairs.exp$gene_id)
  for (j in 1:length(histone_type_list)){
    print(paste("Histone: ", histone_type_list[[j]], sep = ''))
    all_pairs.his = all_pairs.his_list[[j]]
    if (method == "p") {
      all_res_pair = analyze_array(all_pairs.exp, all_pairs.his, option = "p", all_genes=all_genes)
    }
    else if (method == "r") {
      all_res_pair = analyze_array(all_pairs.exp, all_pairs.his, option = "r", n_points=n_points, all_genes=all_genes)
    }
    all_res_list[[j]] = all_res_pair
  }
  return(all_res_list)
}

print("Pearsons-p correlation")

all_res_list.odd_p = analyze_array_list(all_pairs.exp_flt_90_odd, all_pairs.his_list_flt_90_odd, method = "p")
saveRDS(all_res_list.odd_p, "/home/dhthutrang/ENCODE/flank/new_df/all_res_list.odd_p_90_final.RDS")

all_res_list.odd_r = analyze_array_list(all_pairs.exp_flt_90_odd, all_pairs.his_list_flt_90_odd, method = "r")
saveRDS(all_res_list.odd_r, "/home/dhthutrang/ENCODE/flank/new_df/all_res_list.odd_r_90_final.RDS")
