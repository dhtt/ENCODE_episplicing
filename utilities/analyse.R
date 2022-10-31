#LOAD PACKAGES---------
library(fossil)
library(data.table, quietly=TRUE)
library(tidyverse)
library(dplyr, quietly=TRUE)
library(plyr)
# library(rtracklayer)
library(ggplot2)
library(reshape2)
library(NMF)
library(formattable)
# library(clusterProfiler)
# library(org.Hs.eg.db)
# library(DOSE)
library(viridis)
library(ggplot2)
library(ggraph)
library(ggpubr) #plot1 plot2 are in analyze_flank.R 
library(VennDiagram)
library(RColorBrewer)

library(doMC)
doMC::registerDoMC(cores = 8)
#----- GLOBAL VARIABLES -----
setwd("/home/dhthutrang/ENCODE/flank")
workdir = "/home/dhthutrang/ENCODE/flank"
dataset_paths = paste(workdir, c("050722", "110722_0", "110722_1", "110722_2", "190722_0", "190722_1"), sep="/")
dataset_names = c("Abs_Pearson_both.stat", "Abs_Spearman_both.stat", "True_Pearson_both.stat", "True_Spearman_both.stat", "True_Pearson_both.log", "True_Spearman_both.log")
names(dataset_paths) = dataset_names
official_name = c("Adipose tissue", "Aorta", "CD4-positive alpha beta T cell", 
                  "CD8 positive alpha beta T cell", "Ectodermal cell", "Endodermal cell", 
                  "Esophagus", "H1 cell", "Mesenchymal stem cell", "Mesendoderm", "Mesodermal cell", 
                  "Neuronal stem cell", "Pancreas", "Psoas muscle", "Sigmoid colon", 
                  "Small intestine", "Spleen", "Stomach", "Trophoblast")
histone_type_list = c("H3K27ac", "H3K27me3", "H3K36me3", "H3K4me3", "H3K9me3")
subscript = "90_manorm"
epigenomes_annot = read.csv("/home/dhthutrang/ENCODE/utilities/epi_info.csv", 
                            header = TRUE, sep=";", row.names = 1)
#----- ANALYSE METHODS -----
get_p_adj <- function(pearcor_p){
  pearcor_p_adj = vector("list")
  for (i in 1:length(pearcor_p)){
    gene_id = pearcor_p[[i]]$gene_id
    res_list = pearcor_p[[i]]
    res_list_p = apply(res_list[, 2:ncol(res_list)], 1, p.adjust, method = "fdr")
    res_list_p = as.data.frame(cbind(gene_id, t(res_list_p)))
    pearcor_p_adj[[i]] = res_list_p
  }
  names(pearcor_p_adj) = histone_type_list
  return(pearcor_p_adj)
}
filter_by_p <- function(r_df, p_df, option="r") {
  pear_df_filtered_p = list()
  pear_df_filtered_r = list()
  for (i in 1:length(r_df)){
    r_table = as.matrix(r_df[[i]][, 2:ncol(r_df[[i]])])
    p_table = as.matrix(p_df[[i]][, 2:ncol(p_df[[i]])])
    r_table = apply(r_table, 2, as.numeric)
    p_table = apply(p_table, 2, as.numeric)
    if (option == "r"){
      r_table[p_table > 0.05] = NA
      r_table = as.data.frame(r_table)
      r_table = cbind(r_df[[i]][, 1:1], r_table)
      colnames(r_table)[1] = "gene_id"
      pear_df_filtered_r[[i]] = r_table
    }
    else if (option=="p") {
      p_table[abs(r_table) < 0.5] = NA
      p_table = as.data.frame(p_table)
      p_table = cbind(p_df[[i]][, 1:1], p_table)
      colnames(p_table)[1] = "gene_id"
      pear_df_filtered_p[[i]] = p_table
    }
  }
  if (option=="r"){
    names(pear_df_filtered_r) = histone_type_list
    return(pear_df_filtered_r)
  }
  else {
    names(pear_df_filtered_p) = histone_type_list
    return(pear_df_filtered_p)
  }
}
prep_cor_df = function(df_list){
  names(df_list) = histone_type_list
  all_r_his = lapply(histone_type_list, function(x) {
    df = df_list[[x]]
    df = unlist(as.list(df[, 2:length(df)]))
    df = as.data.frame(as.numeric(df[!is.na(df)]))
    df$his = x
    return(df)
  })
  names(all_r_his) = histone_type_list
  all_r_his = do.call(rbind, all_r_his)
  colnames(all_r_his) = c("val", "his") 
  return(all_r_his)
}
get_res_list_sig <- function(res_list, method, r_sig=0.5, p_sig=0.05, exp_sig=2, his_sig=1, compare_threshold = "absolute"){
  res_list_sig = vector("list", length(res_list))
  for (i in 1:length(res_list)) {
    all_res = res_list[[i]]
    
    all_res_sig = vector("list", ncol(all_res)-1)
    for (j in 2:ncol(all_res)){
      all_res.col = as.numeric(all_res[[j]])
      if (method == "pearcor_r"){
        if (compare_threshold == "greater"){
            all_res_sig[[j-1]] = all_res[all_res.col >= r_sig  & !is.na(all_res.col), "gene_id"]
        }
        else if (compare_threshold == "lesser"){
            all_res_sig[[j-1]] = all_res[all_res.col <= r_sig  & !is.na(all_res.col), "gene_id"]
        }
        else if (compare_threshold == "absolute") {
            all_res_sig[[j-1]] = all_res[abs(all_res.col) >= r_sig  & !is.na(all_res.col), "gene_id"]
        }
      }
      else if (method == "pearcor_p" ){
        all_res_sig[[j-1]] = all_res[all_res.col <= p_sig & !is.na(all_res.col), "gene_id"]
      }
      else if (method == "exp" ){
        all_res_sig[[j-1]] = all_res[abs(all_res.col) >= exp_sig & !is.na(all_res.col), "gene_id"]
      }
      else if (method == "his" ){
        all_res_sig[[j-1]] = all_res[abs(all_res.col) >= his_sig & !is.na(all_res.col), "gene_id"]
      }
    }
    names(all_res_sig) = colnames(all_res)[2:length(colnames(all_res))]
    res_list_sig[[i]] = all_res_sig
  }
  return(res_list_sig)
}
get_tissue_spec_array <- function(all_pearcor_p_sig){
  all_tissue_index = vector("list", length(all_pearcor_p_sig))
  for (i in 1:length(all_pearcor_p_sig)){
    all_tissue_name = names(all_pearcor_p_sig[[i]])
    tissue_name = as.vector(sapply(all_tissue_name, function(x) strsplit(x, split = "_")[[1]]))
    tissue_name = unique(tissue_name)
    tissue_index = lapply(tissue_name, function(x) grep(x, all_tissue_name))
    names(tissue_index) = tissue_name
    all_tissue_index[[i]] = tissue_index
  }
  names(all_tissue_index) = histone_type_list
  return(all_tissue_index)
}
get_genes_for_tissue <- function(res_list, tissue_type_list_){
  all_genes_joined = vector("list") 
  for (i in 1:length(res_list)){ 
    res = res_list[[i]]
    tissue_list = tissue_type_list_[[i]]
    all_tissues = lapply(tissue_list, function(x) return(Reduce(union, res[x])))
    names(all_tissues) = names(tissue_list)
    all_genes_joined[[i]] = all_tissues
  }
  names(all_genes_joined) = names(tissue_type_list_)
  return(all_genes_joined)
}
get_genes_for_tissue_filtered <- function(all_genes_joined, TSI_genes){
  all_genes_joined_filter = vector('list')
  for (i in 1:length(all_genes_joined)){
    print(paste("New histone type", histone_type_list[i], sep=' '))
    all_genes_joined_filter[[i]] = vector('list')
    for (j in 1: length(all_genes_joined[[i]])){
      gene_set = all_genes_joined[[i]][[j]]
      print(paste(length(gene_set[gene_set %in% TSI_genes]), length(gene_set)))
      gene_set = gene_set[gene_set %in% TSI_genes]
      print('check')
      all_genes_joined_filter[[i]][[j]] = gene_set
    }
    names(all_genes_joined_filter[[i]]) = names(all_genes_joined[[i]])
  } #List of 6 histone, each has FILTERED sig genes for 25 tissues
  names(all_genes_joined_filter) = names(all_genes_joined)
  return(all_genes_joined_filter)
}
get_all_len <- function(all_genes_joined_list){
  all_len_before = lapply(all_genes_joined_list, function(x) lapply(x, function(y) length(y)))
  all_len_before = lapply(all_len_before, function(x) as.data.frame(cbind(names(x), do.call(rbind, x))))
  all_len_before = all_len_before %>% purrr::reduce(full_join, by = "V1")
  colnames(all_len_before) = c('Tissue', histone_type_list)
  all_len_before$Tissue = official_name
  return(all_len_before)
}

get_overlap_set <- function(list_genes_joined_padj){
  overlap_gene_by_tissue = vector("list")
  for (i in 1:length(names(list_genes_joined_padj$H3K27ac))){
    e_id = names(list_genes_joined_padj$H3K27ac)[i]
    all_res_sig = vector("list")
    for (j in 1:length(list_genes_joined_padj)){
      e_idx = match(e_id, names(list_genes_joined_padj[[j]]))
      print( paste(j, i, e_idx, e_id, sep = ", "))
      all_res_sig[[j]] = list_genes_joined_padj[[j]][[e_idx]]
      }
    gene_set = Reduce(union, all_res_sig)
    overlap_gene_by_tissue[[i]] = gene_set
  }
  return(overlap_gene_by_tissue)
}
compute_jaccard <- function(gene_lists){
  distances = sapply(gene_lists, function(x) sapply(gene_lists, function(y) length(intersect(x, y))))
  distances_union = sapply(gene_lists, function(x) sapply(gene_lists, function(y) length(union(x, y))))
  diag(distances) = NA
  distances = distances/distances_union 
  return(distances)
}
get_sim_matrix <- function(genes_joined_list, name){
  all_sims = vector("list")
  for (i in 1:length(genes_joined_list)){
    print(histone_type_list[i])
    genes_joined_list_his = genes_joined_list[[i]]
    distances = compute_jaccard(genes_joined_list_his) #Jaccard
    
    rownames(distances) = sapply(names(genes_joined_list_his), function(x) epigenomes_annot$Official.name[match(x, rownames(epigenomes_annot))])
    colnames(distances) = rownames(distances)
    
    all_sims[[i]] = distances
  }
  names(all_sims) = name
  print(paste('Highest index: ', max(sapply(all_sims, function(x) max(x, na.rm = T)))))
  print(paste('Lowest index: ', min(sapply(all_sims, function(x) min(x, na.rm = T)))))
  return(all_sims)
}

#----- RUN PIPELINE -----
run_pipeline <- function(dataset_path_, dataset_name_, r_sig=0.5){
  info = list()
  results = list()
  absolute_values = strsplit(dataset_name_, "_")[[1]][1]
  correlation_type = strsplit(dataset_name_, "_")[[1]][2]
  val_type = strsplit(dataset_name_, "_")[[1]][3]
  
  info['absolute_values'] = absolute_values
  info['correlation_type'] = correlation_type
  info['val_type'] = val_type

  print(paste(
    "Currently processing: ", dataset_path_, " with options: Absolute values - ", absolute_values, 
    " Correlation type - ", correlation_type,
    "Value type - ", val_type,
    sep="")
  )


  # Define paths
  print("..... Define paths")
  dataset_path = dataset_path_
  res_p_value_path = paste(paste(dataset_path, "all_res_list.pearcor_p_", sep="/"), subscript,".RDS", sep="")
  res_r_value_path = paste(paste(dataset_path, "all_res_list.pearcor_r_", sep="/"), subscript,".RDS", sep="")
  exp_path = paste(dataset_path, "all_pairs.exp_flt_90.RDS", sep="/")
  his_path = paste(dataset_path, "all_pairs.his_list_flt_90_manorm.RDS", sep="/")

  info['dataset_path'] = dataset_path
  info['res_p_value_path'] = res_p_value_path
  info['res_r_value_path'] = res_r_value_path
  info['exp_path'] = exp_path
  info['his_path'] = his_path


  # Correct R
  print("..... Correct R")
  all_pearcor_p = readRDS(res_p_value_path)
  all_pearcor_padj = get_p_adj(all_pearcor_p)

  results['all_pearcor_padj'] = list(all_pearcor_padj)


  #Filter R values by p-values
  print("..... Filter R values by p-values")
  all_pearcor_r = readRDS(res_r_value_path)
  all_pearcor_r_filtered = filter_by_p(all_pearcor_r, all_pearcor_padj, option = "r")

  results['all_pearcor_r_filtered'] = list(all_pearcor_r_filtered)


  #Get data as controls
  print("..... Get data as controls")
  exp_df = readRDS(exp_path)
  all_pearcor_r_filtered_exp = exp_df %>%
    dplyr::select(-exon_id) %>%
    dplyr::group_by(gene_id) %>%
    dplyr::mutate_all(abs) %>%
    dplyr::summarise_all(sum) %>%
    as.data.frame()
  all_pearcor_r_filtered_exp = list(all_pearcor_r_filtered_exp)

  results['all_pearcor_r_filtered_exp'] = list(all_pearcor_r_filtered_exp)


  his_df = readRDS(his_path)
  all_pearcor_r_filtered_his = lapply(his_df, function(x) 
    return(x %>% dplyr::select(-exon_id) %>%
            dplyr::group_by(gene_id) %>%
            dplyr::mutate_all(abs) %>%
            dplyr::summarise_all(sum) %>%
            as.data.frame()))
  names(all_pearcor_r_filtered_his) = histone_type_list

  results['all_pearcor_r_filtered_his'] = list(all_pearcor_r_filtered_his)


  # Get cumulative plot for R-values
  print("..... Get cumulative plot for R-values")
  all_pearcor_r_cor = prep_cor_df(all_pearcor_r)
  all_pearcor_r_cor$type = "padj > 0.05"
  all_pearcor_r_filtered_cor = prep_cor_df(all_pearcor_r_filtered)
  all_pearcor_r_filtered_cor$type = "padj <= 0.05"
  all_pearcor_r_compare = rbind(all_pearcor_r_cor, all_pearcor_r_filtered_cor)

  results['all_pearcor_r_compare'] = list(all_pearcor_r_compare)

  # Summarize significant genes for each tissue pair
  print("..... Summarize significant genes for each tissue pair")
  all_pearcor_padj_sig = get_res_list_sig(all_pearcor_r_filtered, "pearcor_r", r_sig=r_sig)
  all_pearcor_padj_sig_exp = get_res_list_sig(all_pearcor_r_filtered_exp, "exp", exp_sig=2)
  all_pearcor_padj_sig_his = get_res_list_sig(all_pearcor_r_filtered_his, "his", his_sig=1)

  results['all_pearcor_padj_sig'] = list(all_pearcor_padj_sig)
  results['all_pearcor_padj_sig_exp'] = list(all_pearcor_padj_sig_exp)
  results['all_pearcor_padj_sig_his'] = list(all_pearcor_padj_sig_his)


  # Get significant genes for each tissue
  print("..... Get significant genes for each tissue")
  tissue_type_list = get_tissue_spec_array(all_pearcor_padj_sig)
  tissue_type_list_exp = tissue_type_list[1:1]
  names(tissue_type_list_exp) = "exp"

  all_genes_joined_padj = get_genes_for_tissue(all_pearcor_padj_sig, tissue_type_list)
  all_genes_joined_padj_exp = get_genes_for_tissue(all_pearcor_padj_sig_exp, tissue_type_list_exp)
  all_genes_joined_padj_his = get_genes_for_tissue(all_pearcor_padj_sig_his, tissue_type_list)

  results['all_genes_joined_padj'] = list(all_genes_joined_padj)
  results['all_genes_joined_padj_exp'] = list(all_genes_joined_padj_exp)
  results['all_genes_joined_padj_his'] = list(all_genes_joined_padj_his)


  all_sig_genes = unique(unlist(unlist(all_genes_joined_padj)))
  print(paste("..... ..... Number of total significant genes: ", length(all_sig_genes), sep=""))

  results['all_sig_genes'] = list(all_sig_genes)


  no_sig_gene_by_tissue = get_all_len(all_genes_joined_padj)
  no_sig_gene_by_tissue$`Total epispliced genes (with overlaps)` = rowSums(apply(no_sig_gene_by_tissue[,2:6], 2, as.numeric), na.rm = TRUE)
  no_sig_gene_by_tissue$`Total epispliced genes (without overlaps)` = sapply(get_overlap_set(all_genes_joined_padj), length)
  print(head(no_sig_gene_by_tissue))

  results['no_sig_gene_by_tissue'] = list(no_sig_gene_by_tissue)


  # Get heatmap for distance between tissue pairs
  sim_mat = get_sim_matrix(all_genes_joined_padj, histone_type_list)
  sim_mat_exp = get_sim_matrix(all_genes_joined_padj_exp, "exp")
  sim_mat_his = get_sim_matrix(all_genes_joined_padj_his, histone_type_list)
  sim_mat_exphis = c(sim_mat_exp, sim_mat_his, sim_mat_exp)

  results['sim_mat'] = list(sim_mat)
  results['sim_mat_exphis'] = list(sim_mat_exphis)


  # If there is pos and neg
  if (absolute_values == "True"){
    all_pearcor_padj_sig_pos = get_res_list_sig(all_pearcor_r_filtered, "pearcor_r", r_sig=0, compare_threshold = "greater")
    all_pearcor_padj_sig_neg = get_res_list_sig(all_pearcor_r_filtered, "pearcor_r", r_sig=0, compare_threshold = "lesser")
    
    results['all_pearcor_padj_sig_pos'] = list(all_pearcor_padj_sig_pos)
    results['all_pearcor_padj_sig_neg'] = list(all_pearcor_padj_sig_neg)


    all_genes_joined_padj_pos = get_genes_for_tissue(all_pearcor_padj_sig_pos, tissue_type_list)
    all_genes_joined_padj_neg = get_genes_for_tissue(all_pearcor_padj_sig_neg, tissue_type_list)

    results['all_genes_joined_padj_pos'] = list(all_genes_joined_padj_pos)
    results['all_genes_joined_padj_neg'] = list(all_genes_joined_padj_neg)


    sim_mat_pos = get_sim_matrix(all_genes_joined_padj_pos, histone_type_list)
    sim_mat_neg = get_sim_matrix(all_genes_joined_padj_neg, histone_type_list)

    results['sim_mat_pos'] = list(sim_mat_pos)
    results['sim_mat_neg'] = list(sim_mat_neg)
  }
  else {
    results['all_pearcor_padj_sig_pos'] = NULL
    results['all_pearcor_padj_sig_neg'] = NULL
    results['all_genes_joined_padj_pos'] = NULL
    results['all_genes_joined_padj_neg'] = NULL
    results['sim_mat_pos'] = NULL
    results['sim_mat_neg'] = NULL
  }

  outputs = list(info, results)
  names(outputs) = c("info", "results")
  return(outputs)
}
# results_set_no.R.threshold <- foreach( i=1:length(dataset_paths), .combine='list', .multicombine = TRUE) %dopar% { 
#   output <- run_pipeline(dataset_paths[i], names(dataset_paths)[i], r_sig=0)
# }
# saveRDS(results_set_no.R.threshold, paste(workdir, "results_set_no.R.threshold.RDS", sep='/'))

# results_set <- foreach( i=1:length(dataset_paths), .combine='list', .multicombine = TRUE) %dopar% { 
#   output <- run_pipeline(dataset_paths[i], names(dataset_paths)[i], r_sig=0.5)
# }
# saveRDS(results_set, paste(workdir, "results_set.RDS", sep='/'))
results_set = readRDS(paste(workdir, "results_set.RDS", sep='/'))
print('FIN')
print(names(results_set[[2]][['results']]))

#----- GENERATE PLOTS & TABLES -----
plot_ecdf <- function(result_){
  dataset_path = result_[['info']][['dataset_path']]
  all_pearcor_r_compare = as.data.frame(result_[['results']][['all_pearcor_r_compare']])

  cumulative_r_ecdf = ggplot(all_pearcor_r_compare, 
                           aes(val, colour = his, linetype=type)) + 
    stat_ecdf(aes(colour = his)) +
    xlab("Pearson correlation R values") + ylab("Probability") +
    scale_linetype_manual(values=c("dotted", "solid"), name="", 
                          labels=c(expression(p[adj]<=0.05), expression(p[adj]<=1.0) )) +
    scale_color_discrete(name = "Histone type") +
    theme_minimal() +
    theme(axis.title = element_text(face="bold", size=13), 
          aspect.ratio = 1,
          panel.border = element_rect(colour = "grey", fill=NA, size=0.5),
          legend.position = 'bottom',
          legend.box="vertical", 
          legend.margin=margin()) +
    guides(fill=guide_legend(nrow =3, byrow=t))

  tiff(paste(dataset_path, 'res', paste("cumulative_r_ecdf_", subscript,".tiff", sep=''), sep='/'), width = 8, height = 5, units = 'in', res=300)
  cumulative_r_ecdf
  dev.off()
  print(paste("Plotted ECDF for: ", dataset_path, sep=""))
}
# lapply(results_set, plot_ecdf)

export_no_sig_gene_by_tissue <- function(result_){
  dataset_path = result_[['info']][['dataset_path']]
  no_sig_gene_by_tissue = as.data.frame(result_[['results']][['no_sig_gene_by_tissue']])
  print(no_sig_gene_by_tissue)

  fwrite(
    no_sig_gene_by_tissue,
    paste(dataset_path, 'res', paste("no_sig_gene_by_tissue", subscript,".csv", sep=''), sep='/'), 
    quote=FALSE, row.names = FALSE, col.names = TRUE, sep='\t')
}
# lapply(results_set, export_no_sig_gene_by_tissue)

plot_heatmap <- function(result_){
  dataset_path = result_[['info']][['dataset_path']]
  absolute_values = result_[['info']][['absolute_values']]
  sim_mat = result_[['results']][['sim_mat']]
  sim_mat_exphis = result_[['results']][['sim_mat_exphis']]
  # max_distance = ceiling(max(unlist(lapply(sim_mat, function(x) max(x, na.rm = T)))) * 100)/100
  # min_distance = floor(min(unlist(lapply(sim_mat, function(x) min(x, na.rm = T)))) * 100)/100
  max_distance = 0.5
  min_distance = 0
  print(paste(max_distance, min_distance))

  epigenomes_colors = list(
    list(c("lightgrey", "hotpink", "black"),
        c("plum", "lightblue", "gold", "darkseagreen"),
        c("lightgrey", "hotpink", "black", "MediumPurple"),
        c("plum", "gold")), #1, 2, 5, 6
    list(c("lightgrey", "hotpink", "black"),
        c("plum", "lightblue", "gold", "darkseagreen"),
        c("hotpink", "black", "MediumPurple"),
        c("plum", "gold")) #3, 4
  ) 
  breaks = seq(0, max_distance, (max_distance-min_distance)/101)
  color = '-RdYlBu2:101'

  tiff(
      paste(dataset_path, 'res', paste("heatmap_", subscript,".tiff", sep=''), sep='/'), 
      width = 16, height = 14, units="in", res=200
      )
  par(mfrow = c(2, 3), mar = c(1,1,1,1)*10)
  for (i in 1:(length(sim_mat)+1)){
    if (i %in% c(1,2,3,4,5)){
      aheatmap(sim_mat[[i]], Rowv = FALSE, Colv = TRUE, scale="none", 
              main=paste(histone_type_list[[i]]), 
              annCol = epigenomes_annot[names(rownames(sim_mat[[i]])), c("Potency", "Type", "Origin", "Life.stage")],
              # annRow = epigenomes_annot[names(rownames(sim_mat[[i]])), c("Potency", "Type", "Origin", "Life.stage")],
              annColors = epigenomes_colors[[1]], breaks = breaks, color = color, 
              distfun = 'maximum', hclustfun = "average", 
              reorderfun = function(d, w) reorder(d, w, agglo.FUN = mean), #https://www.rdocumentation.org/packages/gplots/versions/3.1.1/topics/heatmap.2
              # reorderfun = "stats::reorder(dd, 10:1, mean)",
              legend=F, annLegend=F, fontsize=14, cexRow=0, cexCol=1, treeheight=30)
      } 
    else {
      aheatmap(sim_mat[[i-1]], Rowv = FALSE, Colv = TRUE, scale="none",
               main="",
               annCol = epigenomes_annot[names(rownames(sim_mat[[i-1]])), c("Potency", "Type", "Origin", "Life.stage")],
               # annRow = epigenomes_annot[names(rownames(all_sims[[i]])), c("Potency", "Type", "Origin", "Life.stage")],
               annColors = epigenomes_colors[[1]], breaks = breaks, color = color,
               legend=T, annLegend=T, fontsize=14, cexRow=0, cexCol = 1)
    }
  }
  dev.off()


  tiff(
      paste(dataset_path, 'res', paste("heatmap_", subscript,"_exp_his.tiff", sep=''), sep='/'),
      width = 16, height = 14, units="in", res=200
      )
  par(mfrow = c(2, 3), mar = c(1,1,1,1)*10)
  for (i in 1:length(sim_mat_exphis)){
    if (i %in% c(1,2,3,4,5,6)){
      aheatmap(sim_mat_exphis[[i]], Rowv = FALSE, Colv = TRUE, scale="none", 
              main=paste(c('DEU', paste("DHM", histone_type_list, sep='-'), ' ')[i]), 
              annCol = epigenomes_annot[names(rownames(sim_mat_exphis[[i]])), c("Potency", "Type", "Origin", "Life.stage")],
              # annRow = epigenomes_annot[names(rownames(sim_mat[[i]])), c("Potency", "Type", "Origin", "Life.stage")],
              annColors = epigenomes_colors[[1]], color = color, 
              distfun = 'maximum', hclustfun = "average", breaks = seq(0, 1, 1/101),
              reorderfun = function(d, w) reorder(d, w, agglo.FUN = mean), #https://www.rdocumentation.org/packages/gplots/versions/3.1.1/topics/heatmap.2
              # reorderfun = "stats::reorder(dd, 10:1, mean)",
              legend=F, annLegend=F, fontsize=14, cexRow=0, cexCol=1, treeheight=30)
    } 
    else {
      # aheatmap(sim_mat_exphis[[i]], Rowv = FALSE, Colv = TRUE, scale="none", 
      #          main=paste(c('DEU', histone_type_list, ' ')[[i]]), 
      #          annCol = epigenomes_annot[names(rownames(sim_mat_exphis[[i]])), c("Potency", "Type", "Origin", "Life.stage")],
      #          # annRow = epigenomes_annot[names(rownames(sim_mat[[i]])), c("Potency", "Type", "Origin", "Life.stage")],
      #          annColors = epigenomes_colors[[1]], color = color, 
      #          distfun = 'maximum', hclustfun = "complete", 
      #          reorderfun = function(d, w) reorder(d, w, agglo.FUN = mean), #https://www.rdocumentation.org/packages/gplots/versions/3.1.1/topics/heatmap.2
      #          # reorderfun = "stats::reorder(dd, 10:1, mean)",
      #          legend=T, annLegend=T, fontsize=14, cexRow=0, cexCol=1, treeheight=30)
    }
  }
  dev.off()

  if (absolute_values == "True"){
    for (sign_type in c('pos', 'neg')){
      sim_mat_signed = result_[['results']][[paste("sim_mat_", sign_type, sep = '')]]

      tiff(
        paste(dataset_path, 'res', paste("heatmap_", subscript, "_", sign_type, ".tiff", sep=''), sep='/'), 
        width = 16, height = 14, units="in", res=200
      )
      par(mfrow = c(2, 3), mar = c(1,1,1,1)*10)
      for (i in 1:(length(sim_mat_signed)+1)){
        if (i %in% c(1,2,3,4,5)){
          aheatmap(sim_mat_signed[[i]], Rowv = FALSE, Colv = TRUE, scale="none", 
                  main=paste(histone_type_list[[i]]), 
                  annCol = epigenomes_annot[names(rownames(sim_mat_signed[[i]])), c("Potency", "Type", "Origin", "Life.stage")],
                  # annRow = epigenomes_annot[names(rownames(sim_mat[[i]])), c("Potency", "Type", "Origin", "Life.stage")],
                  annColors = epigenomes_colors[[1]], breaks = breaks, color = color, 
                  distfun = 'maximum', hclustfun = "average", 
                  reorderfun = function(d, w) reorder(d, w, agglo.FUN = mean), #https://www.rdocumentation.org/packages/gplots/versions/3.1.1/topics/heatmap.2
                  # reorderfun = "stats::reorder(dd, 10:1, mean)",
                  legend=F, annLegend=F, fontsize=14, cexRow=0, cexCol=1, treeheight=30)
          } 
        else {
          aheatmap(sim_mat_signed[[i-1]], Rowv = FALSE, Colv = TRUE, scale="none",
                  main="",
                  annCol = epigenomes_annot[names(rownames(sim_mat[[i-1]])), c("Potency", "Type", "Origin", "Life.stage")],
                  # annRow = epigenomes_annot[names(rownames(all_sims[[i]])), c("Potency", "Type", "Origin", "Life.stage")],
                  annColors = epigenomes_colors[[1]], breaks = breaks, color = color,
                  legend=T, annLegend=T, fontsize=14, cexRow=0, cexCol = 1)
        }
      }
      dev.off()
    }
  }
}
# lapply(results_set, plot_heatmap)

get_name_list <- function(results_set_){
  name_list_ = lapply(results_set_, function(x) paste(x[['info']][['absolute_values']], x[['info']][['correlation_type']], x[['info']][['val_type']], sep='_'))
  return(name_list_)
}
compare_all_sig_genes <- function(results_set_){
  all_sig_genes_list = lapply(results_set_, function(x) x[['results']][['all_sig_genes']])
  name_list = get_name_list(results_set_) 
  names(all_sig_genes_list) = name_list

  distances = as.data.frame(compute_jaccard(all_sig_genes_list), row.names = name_list)
  colnames(distances) = name_list
  print(distances)
  return(all_sig_genes_list)
}
# all_sig_genes_overlap = compare_all_sig_genes(results_set)
# gene_sets = all_sig_genes_overlap[c(1,2,5,6)]
# lapply(gene_sets, length)
# length(intersect(gene_sets[[3]], gene_sets[[4]]))
# names(gene_sets) = c("Absolute - Pearson", "Absolute - Spearman", "True - Pearson", "True - Spearman")
# myCol <- brewer.pal(length(gene_sets), "Pastel2")
# venn.diagram(
#         x = gene_sets,
#         category.names = names(gene_sets),
#         filename = paste(workdir, 'all_sig_gene_overlaps.png', sep='/'),
#         output=TRUE,
#         imagetype="png" ,
#         height = 5 , 
#         width = 5 , 
#         units = "in",
#         resolution = 300,
#         cat.cex = 0.8,
#         cat.fontface = "bold",
#         cat.default.pos = "outer",
#         lwd = 2,
#         lty = 'blank',
#         fill = myCol,
#         cat.fontfamily = "sans",
#         cex = 0.8,
#         fontfamily = "sans",
#         cat.pos = c(0, 0, 0, 0),
#         print.mode = "percent", 
#         sigdigs = 2
# )

compare_cors <- function(results_set_, cor_method){
  all_cors_list = do.call("cbind", lapply(results_set, function(x) melt(readRDS(x[['info']][['res_r_value_path']]))$value))
  # all_cors_list = do.call("cbind", lapply(results_set_, function(x) melt(as.data.frame(x[['results']][['all_pearcor_r_filtered']][[1]]))$value))
  all_cors_list[is.na(all_cors_list)] = 0
  all_cors_list = all_cors_list[apply(all_cors_list, 1, function(x) !all(x==0)), ]
  name_list = get_name_list(results_set_) 

  cor_mat = cor(all_cors_list, method = cor_method)
  colnames(cor_mat) = name_list
  rownames(cor_mat) = name_list
  cor_mat_abs = cor(abs(all_cors_list), method = cor_method)
  colnames(cor_mat_abs) = name_list
  rownames(cor_mat_abs) = name_list
  print(cor_mat)
  print(cor_mat_abs)
  return(list(cor_mat, cor_mat_abs))
}
# cors_cor = compare_cors(results_set, "spearman")

#----- EXAMPLE GENES -----
get_plot_dfs = function(result_, gene, tissue_pair, selected_his_type){
  all_pairs.exp_flt_90 = readRDS(result_[['info']][['exp_path']])
  all_pairs.his_list_flt_90 = readRDS(result_[['info']][['his_path']])
  return(
    lapply(all_pairs.his_list_flt_90[selected_his_type], function(x)
      df = data.frame(
        DEU=all_pairs.exp_flt_90[[tissue_pair]][all_pairs.exp_flt_90$gene_id == gene], 
        DHM=x[[tissue_pair]][x$gene_id == gene]) 
    )
  )
}
# LMNB1_dfs = get_plot_dfs(results_set[[1]], "LMNB1", "neuronalstemcell_pancreas")
# STIM1_dfs = get_plot_dfs(results_set[[1]], "STIM1", "esophagus_mesenchymalstemcell")
# FGFR2_dfs = lapply(results_set, function(x) get_plot_dfs(x, "FGFR2", "mesenchymalstemcell_sigmoidcolon", seq(1, 5, 1)))
LMNB1_dfs = lapply(results_set, function(x) get_plot_dfs(x, "LMNB1", "neuronalstemcell_pancreas", seq(1, 5, 1)))
cor_df = LMNB1_dfs
gene = "LMNB1"
tissue_type = "Mesenchymal stem cell - Sigmoid colon"
plot_cor = function(dfs, title_prefix, cor_type="together"){
  name_dfs = names(dfs)
  plots = list()
  for (i in 1:length(dfs)){
    df = dfs[[i]]
    plot_ = ggscatter(data=df, x = 'DEU', y = 'DHM',
                           add = "reg.line", 
                           add.params = list(color = "#5D2ED2", fill = "lightgray"),
                           conf.int = TRUE, 
                           font.label = c(8, 'plain', 'red'),
                           title = paste(title_prefix, name_dfs[i], sep='\n')) +
          theme(aspect.ratio = 0.8,
            title = element_text(size=9),
            axis.title = element_text(face = 'bold', size=12))
    if (cor_type == "together"){
      plots[[i]] = plot_ + 
        stat_cor(method = "pearson", col="#EE442F", label.x=min(df$DEU) + 0.5, label.y = max(df$DHM) + 0.8) +
        stat_cor(method = "spearman", col="#007997", label.x=min(df$DEU) + 0.5, label.y = max(df$DHM) + 0.5)
    }
    else {
      plots[[i]] = plot_ + 
        stat_cor(method = cor_type, col="#EE442F", label.x=min(df$DEU) + 0.5, label.y = max(df$DHM) + 0.5)
    }
  }
  names(plots) = name_dfs
  return(plots)
}
cor_plots = list()
cor_plots[[1]] = plot_cor(cor_df[[1]], paste(paste(gene, tissue_type, sep='\n'), paste("Absolute values"), sep="\n"), cor_type="together")
cor_plots[[2]] = plot_cor(cor_df[[3]], paste(paste(gene, tissue_type, sep='\n'), paste("True values"), sep="\n"), cor_type="together") 

plot_list = list(cor_plots[[1]][[3]], cor_plots[[2]][[3]])
tiff("/home/dhthutrang/ENCODE/flank/LMNB1_cor.tiff", width = 7, height = 4, units = "in", res = 300)
ggarrange(plotlist = plot_list, ncol=2, align = 'hv')
dev.off()


DEU_res = fread("/home/dhthutrang/ENCODE/mRNA_seq/dexseqcount/res/neuronalstemcell_pancreas_res.csv")
DEU_res = DEU_res[DEU_res$groupID == "LMNB1",]
DHM_res = fread("/home/dhthutrang/ENCODE/chip_seq/H3K36me3/flank/fl/neuronalstemcell_pancreas.txt.fl.txt")
DHM_res = DHM_res[grep("LMNB1", DHM_res$V9), ]
coef = readRDS("/home/dhthutrang/ENCODE/utilities/plot/neuronalstemcell_pancreas.RData.LMNB1.RDS")
LMNB1_coef = as.data.frame(coef$coeff)
LMNB1_coef$exon = sapply(rownames(LMNB1_coef), function(x) strsplit(x, ':')[[1]][2])
LMNB1_coef$DEU = round(DEU_res$stat*DEU_res[[7]]/abs(DEU_res[[7]]), 1)
LMNB1_coef$DEU_sig = ifelse(DEU_res$padj < 0.05, TRUE, FALSE)
LMNB1_coef$exon_start = seq(1, nrow(LMNB1_coef), 1) - 0.4
LMNB1_coef$exon_end = seq(1, nrow(LMNB1_coef), 1) + 0.4
LMNB1_coef = rbind(LMNB1_coef, NA)
LMNB1_coef$link_to_start_x = c(LMNB1_coef$exon_start[2:nrow(LMNB1_coef)], NA)
LMNB1_coef$DEU_pos = apply(LMNB1_coef[, 1:2], 1, max) + 1

LMNB1_coef_melt = reshape2::melt(LMNB1_coef, value.name = "Exon usage", variable.name = "tissue", 
                        id.vars=c("exon", "exon_start", "exon_end", "link_to_start_x", 
                                  "DEU", "DEU_sig", "DEU_pos"))
LMNB1_coef_melt$link_to_start_y = c(LMNB1_coef_melt$`Exon usage`[2:nrow(LMNB1_coef_melt)], NA)
LMNB1_coef_melt$DHM = unlist(c(DHM_res[1:(nrow(DHM_res)/2), "V10"], NA, DHM_res[(nrow(DHM_res)/2+1):nrow(DHM_res), "V10"], NA))


tiff("/home/dhthutrang/ENCODE/flank/LMNB1_exonusage.tiff", width = 7, height = 4, units = "in", res = 300)
ggplot(LMNB1_coef_melt, aes(color=tissue)) +
  geom_segment(aes(x=exon_start, y=`Exon usage`, xend=exon_end, yend=`Exon usage`), size=1, alpha=0.5) +
  geom_segment(data=LMNB1_coef_melt[LMNB1_coef_melt$sig == TRUE,], aes(x=exon_start, y=`Exon usage`, xend=exon_end, yend=`Exon usage`), size=3) +
  geom_segment(aes(x=exon_end, y=`Exon usage`, xend=link_to_start_x, yend=link_to_start_y), linetype="dotted")  +
  geom_text(aes(label=DEU, x=(exon_start+exon_end)/2, y=DEU_pos), color="black", size=2) +

  geom_segment(aes(x=1, y=0, xend=nrow(LMNB1_coef)-1, yend=0), color="grey") +
  geom_segment(aes(x=exon_start, y=0, xend=exon_end, yend=0), size=5, color="grey") +
  # geom_segment(data=LMNB1_coef_melt[LMNB1_coef_melt$sig == TRUE,], aes(x=exon_start, y=0, xend=exon_end, yend=0), size=5, color="magenta") +

  geom_segment(aes(x=exon_start, y=0, xend=exon_start, yend=DHM)) + 
  scale_color_manual(values=c("#EE442F", "#63ACBE"), name="Tissue", labels=c("Neuronal stem cell", "Pancreas")) +
  scale_x_continuous(breaks=seq(1, nrow(LMNB1_coef)), labels = LMNB1_coef$exon)
dev.off()
# Check for genes
# temp = Reduce(intersect, gene_sets)
# temp = temp[order(temp)]
# temp = results_set[[6]][['results']][['all_pearcor_r_filtered']]
# temp_val = lapply(temp, function(x)  x[x$gene_id == "FGFR2", ])
# temp_id = lapply(temp_val, function(x) {
#   res = x[!is.na(x)]
#   names(res) = colnames(x)[!is.na(x)]
#   return(res)})
# print(temp_id)

# names(results_set[[2]][['results']])
# temp = results_set[[4]][['results']][['all_pearcor_padj']][['H3K27me3']]
# head(temp)
# temp[temp$gene_id == 'FGFR2', 'mesenchymalstemcell_sigmoidcolon']

#----- GO ANNOTATION ANALYSIS -----
# Run on local machine because packages cannot be installed

# universe_genes = readRDS("/home/dhthutrang/ENCODE/mRNA_seq/script/all_DEU_genes.RDS")
# universe_genes = unique(unlist(lapply(universe_genes, function(x) strsplit(x, ';')[[1]][1])))
# print(length(universe_genes))
# strsplit()[[1]]
# combined_genes = unique(unlist(lapply(universe_genes[grep("+", universe_genes, fixed = TRUE)], function(x) strsplit(x, "+", fixed = TRUE)[[1]])))
# single_genes = universe_genes[grep("+", universe_genes, fixed = TRUE, invert = TRUE)]
# print(paste(single_genes, collapse = ', '))
# length(single_genes)




# temp1_his = readRDS(results_set[[3]][['info']][['his_path']])
# temp1_exp = readRDS(results_set[[3]][['info']][['exp_path']])
# head(temp1$H3K36me3)


# temp1_his_ = temp1_his$H3K36me3[temp1_his$H3K36me3$gene_id == "LMNB1", "endodermalcell_neuronalstemcell"]
# temp1_exp_ = temp1_exp[temp1_exp$gene_id == "LMNB1", "endodermalcell_neuronalstemcell"]

# temp2_his_ = temp1_his$H3K36me3[temp1_his$H3K36me3$gene_id == "LMNB1", "neuronalstemcell_pancreas"]
# temp2_exp_ = temp1_exp[temp1_exp$gene_id == "LMNB1", "neuronalstemcell_pancreas"]

# print(cbind(temp1_exp_, temp2_exp_))
# print(cbind(temp1_his_, temp2_his_))
# print(cor.test(temp1_exp_[[1]], temp1_his_[[1]]))
# print(cor.test(temp2_exp_[[1]], temp2_his_[[1]]))


# temp = results_set[[3]][['results']][['all_pearcor_r_filtered']]$H3K36me3
# temp = temp[, c(1, grep("neuronalstemcell", colnames(temp)))]
# temp$score = rowSums(abs(temp[, 2:ncol(temp)]), na.rm=TRUE)
# temp = temp %>% dplyr::select(gene_id, score)
# head(temp)
# temp[which.max(temp$score),]
# fwrite(temp, "/home/dhthutrang/ENCODE/flank/h3k36me3_score.txt", sep = '\t', col.names=FALSE)


