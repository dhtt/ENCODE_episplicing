#LOAD PACKAGES---------
library("fossil")
library(data.table, quietly=TRUE)
library(tidyverse)
library(dplyr, quietly=TRUE)
library(plyr)
library(rtracklayer)
library(ggplot2)
library(reshape2)
library(NMF)
library(formattable)
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
library(viridis)
library(ggplot2)
library(ggraph)
library(ggpubr) #plot1 plot2 are in analyze_flank.R 

#----- GLOBAL VARIABLES------------
setwd("/Users/trangdo/Documents/Episplicing/ENCODE_episplicing/flank") # Local machine
setwd("/home/dhthutrang/ENCODE/flank")
dataset_path = "/home/dhthutrang/ENCODE/flank/110722_2"
workdir = "/home/dhthutrang/ENCODE/flank"
dataset_paths = paste(workdir, c("050722", "110722_0", "110722_2", "110722_2"), sep='/')
dataset_names = c('A_P', "A_S", "T_P", 'T_S')
official_name = c("Adipose tissue", "Aorta", "CD4-positive alpha beta T cell", 
                  "CD8 positive alpha beta T cell", "Ectodermal cell", "Endodermal cell", 
                  "Esophagus", "H1 cell", "Mesenchymal stem cell", "Mesendoderm", "Mesodermal cell", 
                  "Neuronal stem cell", "Pancreas", "Psoas muscle", "Sigmoid colon", 
                  "Small intestine", "Spleen", "Stomach", "Trophoblast")
histone_type_list = c("H3K27ac", "H3K27me3", "H3K36me3", "H3K4me3", "H3K9me3")
subscript = '90_manorm'
epigenomes_annot = read.csv("/Users/trangdo/Documents/Episplicing/ENCODE_episplicing/utilities/epi_info.csv", 
                            header = TRUE, sep=';', row.names = 1) # Local machine
epigenomes_annot = read.csv("/home/dhthutrang/ENCODE/utilities/epi_info.csv", 
                            header = TRUE, sep=';', row.names = 1)
#----- READ CORRELATION RESULTS----------
# Last version before BMC
res_p_value_path = paste("/Users/trangdo/Documents/Episplicing/ENCODE_episplicing/utilities/corrected_cor/all_res_list.pearcor_p_", subscript,".RDS", sep='')
res_r_value_path = paste("/Users/trangdo/Documents/Episplicing/ENCODE_episplicing/utilities/corrected_cor/all_res_list.pearcor_r_", subscript, ".RDS", sep='')
exp_path = "/Users/trangdo/Documents/Episplicing/ENCODE_episplicing/utilities/manorm_90_get_cor_dataset/all_pairs.exp_flt_90.RDS"
his_path = "/Users/trangdo/Documents/Episplicing/ENCODE_episplicing/utilities/manorm_90_get_cor_dataset/all_pairs.his_list_flt_90_manorm.RDS"

# BMC revision
get_paths <- function(dataset_path_){
  res_p_value_path = paste(paste(dataset_path_, "all_res_list.pearcor_p_", sep='/'), subscript,".RDS", sep='')
  res_r_value_path = paste(paste(dataset_path_, "all_res_list.pearcor_r_", sep='/'), subscript,".RDS", sep='')
  exp_path = paste(dataset_path_, "all_pairs.exp_flt_90.RDS", sep='/')
  his_path = paste(dataset_path_, "all_pairs.his_list_flt_90_manorm.RDS", sep='/')
  paths_ = list(res_p_value_path, res_r_value_path, exp_path, his_path)
  names(paths_) = c('res_p_value_path', 'res_r_value_path', 'exp_path', 'his_path')
  return(paths_)
}

paths_ = lapply(dataset_paths, get_paths)
names(paths_) = dataset_names
#----- Get padj ------- 
#' Compute p_adj for all p_val from pairwise correlation for each histone
#' Return: List of DF with p_adj values
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
all_pearcor_p = readRDS(res_p_value_path)
all_pearcor_padj = get_p_adj(all_pearcor_p)
head(all_pearcor_p[[1]])

# saveRDS(all_pearcor_padj, "all_pearcor_padj.RDS")
table(all_pearcor_padj$H3K36me3 <= 0.05)
temp = all_pearcor_p$H3K36me3
names(all_pearcor_p) = histone_type_list
temp1 = all_pearcor_padj$H3K36me3

#---------- [optional] Filter r_val by padj ------- 
#' Retain r_val if p_adj <= 0.05
#' Return: List of r_val where values are retained for p_adj < 0.05 
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
      colnames(r_table)[1] = 'gene_id'
      pear_df_filtered_r[[i]] = r_table
    }
    else if (option=="p") {
      p_table[abs(r_table) < 0.5] = NA
      p_table = as.data.frame(p_table)
      p_table = cbind(p_df[[i]][, 1:1], p_table)
      colnames(p_table)[1] = 'gene_id'
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
all_pearcor_r = readRDS(res_r_value_path)
# all_pearcor_r = filter_by_p(all_pearcor_r, all_pearcor_padj)
all_pearcor_r_filtered = filter_by_p(all_pearcor_r, all_pearcor_padj, option = 'r')
# all_pearcor_p_filtered = filter_by_p(all_pearcor_r, all_pearcor_padj, option = 'p')

#Get data as controls
exp_df = readRDS(exp_path)
head(exp_df)
all_pearcor_r_filtered_exp = exp_df %>%
  dplyr::select(-exon_id) %>%
  dplyr::group_by(gene_id) %>%
  dplyr::mutate_all(abs) %>%
  dplyr::summarise_all(sum) %>%
  as.data.frame()
all_pearcor_r_filtered_exp = list(all_pearcor_r_filtered_exp)

his_df = readRDS(his_path)
head(his_df$H3K27ac[1:10, 1:10])
all_pearcor_r_filtered_his = lapply(his_df, function(x) 
  return(x %>% dplyr::select(-exon_id) %>%
           dplyr::group_by(gene_id) %>%
           dplyr::mutate_all(abs) %>%
           dplyr::summarise_all(sum) %>%
           as.data.frame()))
names(all_pearcor_r_filtered_his) = histone_type_list

#======TESTING======
lapply(all_pearcor_r_filtered, function(x) head(x[, 1:20]))
View(all_pearcor_r_filtered[[1]])
x = all_pearcor_r_filtered$H3K9me3[!is.na(all_pearcor_r_filtered$H3K9me3$adiposetissue_aorta), ]
x = t(all_pearcor_r_filtered[[6]][all_pearcor_r_filtered$H3K27ac$gene_id == "STIM1",])
rownames(x)[!is.na(x)]
rownames(x) = colnames(all_pearcor_r_filtered$H3K27me3)
table(x)
endodermalcell_mesenchymalstemcell

y = all_pearcor_r_filtered$H3K36me3[all_pearcor_r_filtered$H3K36me3$gene_id == "LMNB1",]
y = all_pearcor_r_filtered$H3K4me3[all_pearcor_r_filtered$H3K4me3$gene_id == "LMNB1",]
y = all_pearcor_r_filtered$H3K27ac[all_pearcor_r_filtered$H3K27ac$gene_id == "LMNB1",]
y[grep("neuron", names(y))]
names(y)[!is.na(y)]
colnames(all_pearcor_r_filtered$H3K36me3)
temp = lapply(all_pearcor_r_filtered, function(y){
  z = y[y$gene_id == "LMNB1",] 
  return(names(z)[!is.na(z)])
})

temp = temp[sapply(temp, function(y) length(y) != 1)]
table(Reduce(c, temp))[order(table(Reduce(c, temp)))]
lapply(temp, function(y) table(unlist(strsplit(y, split='_'))))
temp$H3K9me3[grep('neuronal', temp$H3K9me3)]

common_pairs = Reduce(intersect, lapply(all_pearcor_r_filtered, function(x) colnames(x)))
do.call(sum, lapply(all_pearcor_r_filtered, function(x) apply(x[, common_pairs], 1, function(y) length(y[!is.na(y)]))))
x = data.frame(all_assoc = Reduce(sum, lapply(all_pearcor_r_filtered, function(x) 
  apply(x[, common_pairs], 1, function(y) length(y[!is.na(y)])))),
           names = all_pearcor_r_filtered$H3K27ac$gene_id)
x1 = t(do.call(rbind, lapply(all_pearcor_r_filtered, function(x) print(x[x$gene_id == "TCTN3", common_pairs]))))
Reduce(sum, lapply(all_pearcor_r_filtered, function(x) 
  apply(x[, common_pairs], 1, function(y) length(y[!is.na(y)]))))
x = lapply(all_pearcor_r_filtered, function(x) 
  apply(x[, common_pairs], 1, function(y) length(y[!is.na(y)])))
x = as.data.frame(do.call(cbind, x))
x$sum = rowSums(x)
x$gene = all_pearcor_r_filtered$H3K27ac$gene_id
x1 = x$H3K9me3
table(sapply(x1, length))
x1[sapply(x1, length) == 15]
#----ECDF----
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
all_pearcor_r_cor = prep_cor_df(all_pearcor_r)
all_pearcor_r_cor$type = 'padj > 0.05'
all_pearcor_r_filtered_cor = prep_cor_df(all_pearcor_r_filtered)
all_pearcor_r_filtered_cor$type = 'padj <= 0.05'
all_pearcor_r_compare = rbind(all_pearcor_r_cor, all_pearcor_r_filtered_cor)

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

# all_pearcor_r_cor = prep_cor_df(all_pearcor_padj)
# all_pearcor_r_cor$type = '|r| > 0.5'
# all_pearcor_r_filtered_cor = prep_cor_df(all_pearcor_p_filtered)
# all_pearcor_r_filtered_cor$type = '|r| <= 0.05'
# all_pearcor_r_compare = rbind(all_pearcor_r_cor, all_pearcor_r_filtered_cor)

# cumulative_p_ecdf = ggplot(all_pearcor_r_compare, aes(val, colour = his, linetype=type)) + 
#   stat_ecdf(aes(colour = his)) +
#   xlab("Pearson correlation p-adj values") + ylab("Probability") +
#   scale_linetype_manual(values=c("solid", "dotted", "twodash")) +
#   theme_bw() + theme(aspect.ratio=1)
# cumulative_p_ecdf
#----- Get significant results -----
#' Get entries with p_adj < 0.05
get_res_list_sig <- function(res_list, method, r_sig=0.5, p_sig= 0.05, exp_sig=2, his_sig=1, compare_threshold = "absolute"){
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
        else {
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
# all_pearcor_padj_sig = get_res_list_sig(all_pearcor_r_filtered, "pearcor_r", r_sig=0.5, compare_threshold = "greater")
# all_pearcor_padj_sig = get_res_list_sig(all_pearcor_r_filtered, "pearcor_r", r_sig=-0.5, compare_threshold = "lesser")
all_pearcor_padj_sig = get_res_list_sig(all_pearcor_r_filtered, "pearcor_r", r_sig=-0.5)
all_pearcor_padj_sig_exp = get_res_list_sig(all_pearcor_r_filtered_exp, "exp", exp_sig=2)
all_pearcor_padj_sig_his = get_res_list_sig(all_pearcor_r_filtered_his, "his", his_sig=1)

lapply(all_pearcor_padj_sig_exp, function(x) paste(x$neuronalstemcell_spleen, collapse = ','))
all_pearcor_padj_sig[[1]]$stomach_trophoblast
all_pearcor_padj_sig_exp[[1]]$adiposetissue_aorta
all_sig_genes = Reduce(union, Reduce(union, all_pearcor_padj_sig[[3]]))
length(all_sig_genes)
paste(all_sig_genes, collapse = ', ')
temp = as.data.frame(all_sig_genes)

#----- Get significant genes for each tissue-----
get_tissue_spec_array <- function(all_pearcor_p_sig){
  all_tissue_index = vector("list", length(all_pearcor_p_sig))
  for (i in 1:length(all_pearcor_p_sig)){
    all_tissue_name = names(all_pearcor_p_sig[[i]])
    tissue_name = as.vector(sapply(all_tissue_name, function(x) strsplit(x, split = '_')[[1]]))
    tissue_name = unique(tissue_name)
    tissue_index = lapply(tissue_name, function(x) grep(x, all_tissue_name))
    names(tissue_index) = tissue_name
    all_tissue_index[[i]] = tissue_index
  }
  names(all_tissue_index) = histone_type_list
  return(all_tissue_index)
}
tissue_type_list = get_tissue_spec_array(all_pearcor_padj_sig)
tissue_type_list_exp = tissue_type_list[1:1]
names(tissue_type_list_exp) = 'exp'

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
all_genes_joined_padj = get_genes_for_tissue(all_pearcor_padj_sig, tissue_type_list)
all_genes_joined_padj_exp = get_genes_for_tissue(all_pearcor_padj_sig_exp, tissue_type_list_exp)
all_genes_joined_padj_his = get_genes_for_tissue(all_pearcor_padj_sig_his, tissue_type_list)
paste(all_genes_joined_padj$H3K27me3$neuronalstemcell, collapse = ',')
all_sig_genes = unique(unlist(unlist(all_genes_joined_padj)))
# saveRDS(all_sig_genes, paste("all_sig_genes_", subscript, ".RDS", xep=''))

paste(all_sig_genes, collapse = ', ')
lapply(histone_type_list, function(x) length(unique(unlist(all_genes_joined_padj[[x]]))))

#---------- [optional] Filter by TSI ----- 
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
all_genes_joined_padj_filtered = get_genes_for_tissue_filtered(all_genes_joined_padj, TSG_gene)

#----- Summarize results overlap -----
get_all_len <- function(all_genes_joined_list){
  all_len_before = lapply(all_genes_joined_list, function(x) lapply(x, function(y) length(y)))
  all_len_before = lapply(all_len_before, function(x) as.data.frame(cbind(names(x), do.call(rbind, x))))
  all_len_before = all_len_before %>% purrr::reduce(full_join, by = "V1")
  colnames(all_len_before) = c('Tissue', histone_type_list)
  all_len_before$Tissue = official_name
  return(all_len_before)
}
no_sig_gene_by_tissue = get_all_len(all_genes_joined_padj)
no_sig_gene_by_tissue$`Total epispliced genes (with overlaps)` = rowSums(apply(no_sig_gene_by_tissue[,2:6], 2, as.numeric), na.rm = TRUE)

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
no_sig_gene_by_tissue$`Total epispliced genes (without overlaps)` = sapply(get_overlap_set(all_genes_joined_padj), length)

saveRDS(
    no_sig_gene_by_tissue, 
    paste(dataset_path, 'res', paste("no_sig_gene_by_tissue_", subscript,".RDS", sep=''), sep='/'))
fwrite(
    no_sig_gene_by_tissue,
    paste(dataset_path, 'res', paste("no_sig_gene_by_tissue", subscript,".csv", sep=''), sep='/'), 
    quote=FALSE, row.names = FALSE, col.names = TRUE, sep='\t')

# Formattable
tissue_formatter <- formatter("span", 
                              style = x ~ style(
                              width = suffix(x, "px"),
                              font.weight = "bold", 
                              color = ifelse(x == "Total", "black", "gray")))
#!!!!!!!!!!!!!!!!!! Adipose tissue: 120
# no_sig_gene_by_tissue$`Total epispliced genes (with overlaps)`[1] = min(no_sig_gene_by_tissue$`Total epispliced genes (with overlaps)`)
# no_sig_gene_by_tissue$`Total epispliced genes (without overlaps)`[1] = min(no_sig_gene_by_tissue$`Total epispliced genes (without overlaps)`)
formattable(no_sig_gene_by_tissue,
            # align =c("l", "l", "c","c","c","c"), 
            list(
              # `H3K4me1` = color_tile("white", "wheat"),
              `H3K4me3` = color_tile("white", "wheat"),
              `H3K9me3` = color_tile("white", "wheat"),
              `H3K27me3` = color_tile("white", "wheat"),
              `H3K36me3` = color_tile("white", "wheat"), 
              `H3K4me1` = color_tile("white", "wheat"),
              `H3K27ac` = color_tile("white", "wheat"),
              `Total epispliced genes (without overlaps)` = color_tile("white", "lightsalmon"),
              `Total epispliced genes (with overlaps)` = color_tile("white", "lightsalmon")
            )
)

formattable(no_sig_gene_by_tissue_filtered,
            # align =c("l", "l", "c","c","c","c"), 
            list(
              # `H3K4me1` = color_tile("white", "wheat"),
              `H3K4me3` = color_tile("white", "wheat"),
              `H3K9me3` = color_tile("white", "wheat"),
              `H3K27me3` = color_tile("white", "wheat"),
              `H3K36me3` = color_tile("white", "wheat"), 
              `H3K4me1` = color_tile("white", "wheat"),
              `H3K27ac` = color_tile("white", "wheat"),
              `Total epispliced genes (without overlaps)` = color_tile("white", "lightsalmon"),
              `Total epispliced genes (with overlaps)` = color_tile("white", "lightsalmon")
            )
)

#----- Heatmap of tissue-spec genes -----
get_sim_matrix <- function(genes_joined_list, name){
  all_sims = vector("list")
  for (i in 1:length(genes_joined_list)){
    print(histone_type_list[i])
    genes_joined_list_his = genes_joined_list[[i]]
    distances = sapply(genes_joined_list_his, function(x) sapply(genes_joined_list_his, function(y) length(intersect(x, y))))
    distances_union = sapply(genes_joined_list_his, function(x) sapply(genes_joined_list_his, function(y) length(union(x, y))))
    diag(distances) = NA
    distances = distances/distances_union #Jaccard
    
    rownames(distances) = sapply(names(genes_joined_list_his), function(x) epigenomes_annot$Official.name[match(x, rownames(epigenomes_annot))])
    colnames(distances) = rownames(distances)
    
    all_sims[[i]] = distances
  }
  names(all_sims) = name
  print(paste('Highest index: ', max(sapply(all_sims, function(x) max(x, na.rm = T)))))
  print(paste('Lowest index: ', min(sapply(all_sims, function(x) min(x, na.rm = T)))))
  return(all_sims)
}
sim_mat = get_sim_matrix(all_genes_joined_padj, histone_type_list)
sim_mat_exp = get_sim_matrix(all_genes_joined_padj_exp, 'exp')
sim_mat_his = get_sim_matrix(all_genes_joined_padj_his, histone_type_list)
sim_mat_exphis = c(sim_mat_exp, sim_mat_his, sim_mat_exp)

max_distance = ceiling(max(unlist(lapply(sim_mat, function(x) max(x, na.rm = T)))) * 100)/100
min_distance = floor(min(unlist(lapply(sim_mat, function(x) min(x, na.rm = T)))) * 100)/100

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
for (i in 1:length(sim_mat)){
  if (i %in% c(1,2,3,4,5)){
    aheatmap(sim_mat[[i]], Rowv = FALSE, Colv = TRUE, scale="none", 
             main=paste(histone_type_list[[i]]), 
             annCol = epigenomes_annot[names(rownames(sim_mat[[i]])), c("Potency", "Type", "Origin", "Life.stage")],
             # annRow = epigenomes_annot[names(rownames(sim_mat[[i]])), c("Potency", "Type", "Origin", "Life.stage")],
             annColors = epigenomes_colors[[1]], breaks = breaks, color = color, 
             distfun = 'maximum', hclustfun = "complete", 
             reorderfun = function(d, w) reorder(d, w, agglo.FUN = mean), #https://www.rdocumentation.org/packages/gplots/versions/3.1.1/topics/heatmap.2
             # reorderfun = "stats::reorder(dd, 10:1, mean)",
             legend=F, annLegend=F, fontsize=14, cexRow=0, cexCol=1, treeheight=30)
    } 
  else {
    # aheatmap(sim_mat[[i]], Rowv = FALSE, Colv = TRUE, scale="none",
    #          main=paste(histone_type_list[[i]]),
    #          annCol = epigenomes_annot[names(rownames(sim_mat[[i]])), c("Potency", "Type", "Origin", "Life.stage")],
    #          # annRow = epigenomes_annot[names(rownames(all_sims[[i]])), c("Potency", "Type", "Origin", "Life.stage")],
    #          annColors = epigenomes_colors[[1]], breaks = breaks, color = color,
    #          legend=T, annLegend=T, fontsize=14, cexRow=0, cexCol = 1)
  }
}
dev.off()


tiff(
    paste(dataset_path, 'res', paste("heatmap_", subscript,"_exp_his.tiff", sep=''), sep='/'),
    width = 16, height = 14, units="in", res=200)
par(mfrow = c(2, 3), mar = c(1,1,1,1)*10)
for (i in 1:length(sim_mat_exphis)){
  if (i %in% c(1,2,3,4,5,6)){
    aheatmap(sim_mat_exphis[[i]], Rowv = FALSE, Colv = TRUE, scale="none", 
             main=paste(c('DEU', paste("DHM", histone_type_list, sep='-'), ' ')[[i]]), 
             annCol = epigenomes_annot[names(rownames(sim_mat_exphis[[i]])), c("Potency", "Type", "Origin", "Life.stage")],
             # annRow = epigenomes_annot[names(rownames(sim_mat[[i]])), c("Potency", "Type", "Origin", "Life.stage")],
             annColors = epigenomes_colors[[1]], color = color, 
             distfun = 'maximum', hclustfun = "complete", 
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


#----- Rand Index -----
all_distances = vector("list")
for (i in 1:length(sim_mat)){
  all_distances[[i]] = 1- sim_mat[[i]]
}
reorderfun = function(d, w) reorder(d, w)
reorderfun2 = function(d, w) reorder(d, w, mean)
reorderfun3 = function(d, w) reorder(d, 10:1, mean)

d = dist(all_distances[[1]]) #Change idx for each histone
h = hclust(d, "complete")
d_de = as.hclust(reorderfun(as.dendrogram(h), h$order))
d_de2 = as.hclust(reorderfun2(as.dendrogram(h), h$order))
d_de3 = as.hclust(reorderfun3(as.dendrogram(h), h$order))
plot(d_de)
plot(d_de2)
plot(d_de3)
rect.hclust(d_de, k=2)

#---------- H cluster ----
Hs = list( #H2K27ac
  list(
    c(rep(1,11), rep(2,6), 3, 2),
    c(rep(1,9), rep(2,2), rep(3,6), 4, 3), 
    c(rep(1,2), 2, rep(1,3), rep(2,3), rep(1,2), rep(3,4), 1, rep(3,2), 4),
    c(rep(1,11), rep(2,8))
  ), #H3K27me3
  list(
    c(rep(1,10), rep(2,4), 3, rep(2,2)),
    c(rep(1,7), rep(2,2), 1, rep(3,4), 4, rep(3,2)),
    c(rep(1,2), rep(2,2), rep(1,2), 2, rep(1,2), 2, 3, 1, rep(3,4), 4),
    c(rep(1,10), rep(2,7))
  ), #H3K36me3
  list(
    c(rep(1,8), 2, rep(1,2), rep(2,5), 3),
    c(rep(1,8), 3, rep(2,2), rep(3,5), 4),
    c(1, 2, rep(1,3), rep(2,3), rep(1,3), rep(3,2), 4, rep(3,3)),
    c(rep(1,8), 2, c(1,2), rep(2,6))
  ), #H3K4me3
  list(
    c(rep(1,10), rep(2,6), 3),
    c(rep(1,8), rep(2,2), rep(3,6), 4),
    c(rep(1,2), 2, 1, rep(2,2), 1, 2, rep(1,2), 4, 1, rep(3,5)),
    c(rep(1,10), rep(2,7))
  ), #H3K9me3
  list(
    c(rep(1,8), rep(2, 3), rep(1,2), rep(2,3), 3),
    c(rep(1,8), rep(2, 3), rep(3,2), rep(2,3), 4),
    c(rep(1,2), 2, rep(1,2), rep(2,3), rep(3,3), rep(1,3), 4, rep(3,2)),
    c(rep(1,8), rep(2,3), rep(1,2), rep(2,4))
  )
) 

#--------- Initial_clusters ----
Hs = list(#exp
  list(
    c(rep(1,7), rep(2,2), rep(1,3), rep(2,4), 3, 1, 2),
    c(rep(1,7), rep(2,2), rep(3,2), 1, rep(2,4), 4, 1, 2),
    c(1, 2, 1, rep(2,3), 1, 2, 3, rep(2,3), 3, 4, rep(3,3), 1, 3),
    c(rep(1,7), rep(2,2), rep(1,3), rep(2,5), 1, 2)
  ), #H2K27ac
  list(
    c(rep(1,5), 2,1, rep(2,4), 3, 1, rep(2,6)),
    c(rep(1,5), 2,1, rep(2,4), 3, 1, rep(2,4), rep(3,2)), 
    c(1, 2, rep(1,3), 3, 1, rep(3,2), rep(4,2), 1, rep(4,3), 3, rep(4,3)),
    c(rep(1,5), 2, 1, rep(2,4), rep(1,2), rep(2,6))
  ), #H3K27me3
  list(
    c(rep(1,4), 2, 1, 2, 3, rep(2,3), 1, rep(2,5)),
    c(rep(1,4), 2, 1, 3, 4, rep(2,3), 1, rep(2,4), 3),
    c(1, 2, rep(1,2), 3, rep(4,2), 1, 3, rep(4,2), 1, 4, rep(3,2), rep(4,2)),
    c(rep(1,4), 2, 1, 2, 1,rep(2,3), 1, rep(2,5))
  ), #H3K36me3
  list(
    c(rep(1,2), rep(2,3), 1, rep(2,4), 3, rep(1,2), rep(2,2), 1, 2),
    c(rep(1,2), rep(2,3), 1, rep(2,3), 3, 4, rep(1,2), 3, 2, 1, 2),
    c(1, 2, 3, 2, 3, 1, 3, rep(2,3), rep(1,3), 2, 3, 4, 2),
    c(rep(1,2), rep(2,3), 1, rep(2,4), rep(1,3), rep(2,2), 1, 2)
  ), #H3K4me3
  list(
    c(rep(1,3), rep(2,5), 1, 2, 3, 2, rep(1,2), rep(1,3)),
    c(rep(1,3), rep(2,2), 3, 2, 3, 1, 2, 4, 2, rep(1,2), rep(2,3)),
    c(1,2,3,2,4,2,4,2,1,2,1,2,rep(1,2), rep(4,2), 2),
    c(rep(1,3), rep(2,5), 1, 2, 1, 2, rep(1,2), rep(2,3))
  ), #H3K9me3
  list(
    c(1, rep(2,3), rep(1,2), rep(2,4), 1, 3, rep(2,2), rep(1,2), 2),
    c(1, rep(2,3), rep(1,2), rep(2,4), 1, 3, rep(4,2), rep(1,2), 2),
    c(1, rep(1,2), 3, 4, 3, 2, 3, 2, 3, rep(1,2), rep(3,2), rep(1,2), 3),
    c(1, rep(2,3), rep(1,2), rep(2,4), rep(1,2), rep(2,2), rep(1,2), 2)
  )
) 
clusters = list(#exp
  list(
    c(rep(1,3), rep(2,16)),
    c(rep(1,3), rep(2,9), rep(3,7)),
    c(rep(1,3), rep(2,9), rep(3,7)),
    c(rep(1,3), rep(2,16))
  ), #H1
  list(
    c(rep(1,4), rep(2,15)),
    c(rep(1,4), rep(2,13), rep(3,2)),
    c(rep(1,4), rep(2,13), rep(3,2)),
    c(rep(1,4), rep(2,15))
  ), #H2
  list(
    c(rep(1,5), rep(2,12)),
    c(rep(1,5), rep(2,6), rep(3,6)),
    c(rep(1,5), rep(2,6), rep(3,6)),
    c(rep(1,5), rep(2,12))
  ), #H3
  list(
    c(rep(1,9), rep(2,8)),
    c(rep(1,6), rep(2,3), rep(3,8)),
    c(rep(1,6), rep(2,3), rep(3,8)), 
    c(rep(1,9), rep(2,8))
  ), #H4
  list(
    c(rep(1,5), rep(2,12)),
    c(rep(1,5), rep(2,3), rep(3,9)),
    c(rep(1,5), rep(2,3), rep(3,9)),
    c(rep(1,5), rep(2,12))
  ), #H5
  list(
    c(rep(1,10), rep(2,7)),
    c(rep(1,4), rep(2,6), rep(3,7)),
    c(rep(1,4), rep(2,6), rep(3,7)),
    c(rep(1,10), rep(2,7))
  )
)

#--------- Rand adj ----
library(fossil)
get_rand_adj <- function(Hs_, clusters_, name){
  rand.adj_res = list()
  for (i in 1:length(Hs_)){
    H_rand.adj = c()
    for (j in 1:4){
      rand_i = round(adj.rand.index(clusters_[[i]][[j]], Hs_[[i]][[j]]), 3)
      H_rand.adj = c(H_rand.adj, rand_i)
    }
    rand.adj_res[[i]] = H_rand.adj
  }
  names(rand.adj_res) = name
  rand.adj_res = as.data.frame(rand.adj_res)
  rownames(rand.adj_res) = rev(c("Life stage", "Origin", "Type", "Potency"))
  # rand.adj_res = rand.adj_res[c( "Potency", "Type", "Origin", "Life stage"),]
  return(as.data.frame(rand.adj_res))
} 
rand_adj = get_rand_adj(Hs, clusters, histone_type_list)
rand_adj = get_rand_adj(Hs, clusters, c("DEU", histone_type_list))
homogeneity(Hs[[5]][[2]], clusters[[5]][[2]])

homo = matrix(0, ncol = length(Hs), nrow = length(Hs[[1]]))
temp_row = list()
for (i in 1:length(Hs)){
  temp_col = list()
  for (j in 1:length(Hs[[1]])){
    temp_col[j] = homogeneity(Hs[[i]][[j]], clusters[[i]][[j]])
  }
  temp_row[[i]] = temp_col
}

do.call(cbind, lapply(temp_row, unlist))

fwrite(rand_adj, paste("rand_index_", subscript,".csv", sep=''), sep='\t', quote=FALSE, dec = ',', row.names = TRUE, col.names = TRUE)

#----- Gene annotation -----
get_entrez_id <-function(gene_list){
  entrez_gene_list = AnnotationDbi::select(org.Hs.eg.db,
                            keys = gene_list,
                            columns = c("ENTREZID", "SYMBOL"),
                            keytype = "SYMBOL")
  return(entrez_gene_list[[2]])
}
all_genes_joined_padj_his = lapply(all_genes_joined_padj, function(x) unique(unlist(x)))
gene_list_his = lapply(all_genes_joined_padj_his, get_entrez_id)
names(gene_list_his) = histone_type_list
lapply(all_genes_joined_padj_his, length)
paste(all_genes_joined_padj_his$H3K9me3, collapse = ',')

ck_bp_005_gene = compareCluster(geneCluster = gene_list_his, fun = "enrichGO",
                                OrgDb='org.Hs.eg.db', ont = "BP", qvalueCutoff = 0.05,
                                pAdjustMethod = "fdr", readable =TRUE)
ck_bp_001_gene = compareCluster(geneCluster = gene_list_his, fun = "enrichGO",
                                OrgDb='org.Hs.eg.db', ont = "BP", qvalueCutoff = 0.01,
                                pAdjustMethod = "fdr", readable =TRUE)


annot_df = ck_bp_005_gene
head(annot_df)
annot_df@compareClusterResult$Description = str_to_sentence(annot_df@compareClusterResult$Description)
annot_df@compareClusterResult$Categories = 'Cellular and Metabolic\nProcesses'
annot_df@compareClusterResult$Categories[grep('size|divis|devel|genesis|growth|morpho|matur|organiz|guidan|organel|anato|exten', annot_df@compareClusterResult$Description, ignore.case = T)] = "Developmental\nProcesses"
annot_df@compareClusterResult$Categories[grep('signal|transport|import|secret|export|exo|transmi', annot_df@compareClusterResult$Description, ignore.case = T)] = "Cell Signaling"

annot_df@compareClusterResult$Description = gsub('h3-k9', 'H3K9', annot_df@compareClusterResult$Description)
annot_df@compareClusterResult$Description = gsub('Ncrna', 'ncRNA', annot_df@compareClusterResult$Description)
annot_df@compareClusterResult$Description = gsub('Snrna', 'snRNA', annot_df@compareClusterResult$Description)
# tiff(paste("new_fig/final/annot_", subscript,".tiff", sep=''), height = 10, width = 10, res=300, units = 'in')

dotplot = clusterProfiler::dotplot(annot_df, showCategory = 12) +
  facet_wrap(vars(Categories), nrow=2, scales = "free") +
  scale_color_viridis(option = "D", name="FDR-adjusted p-value") +
  scale_size(name='Gene Ratio') +
  theme(
    plot.title = element_text(size = 10, hjust = 0.5, face = "bold"),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    axis.text.x = element_text(colour="black", size = 10, color = 'black', angle = 45, vjust=0.5),
    axis.text.y = element_text(colour="black", size = 12, color = 'black'),
    panel.background = element_rect(fill = 'white'),
    legend.position = c("right"),
    legend.box="vertical", 
    legend.margin=margin(),
    legend.justification = "top",
    strip.background = element_rect(fill="#5D2ED2"),
    strip.text = element_text(color = 'white', face='bold', size=12)
    )
tiff(paste("/Users/trangdo/Google\ Drive/manuscript_figures/figures/manorm_90/annot_", subscript,"_3.tiff", sep=''), 
     height = 11, width = 14, res=400, units = 'in')
# tiff(paste("/Users/trangdo/Google\ Drive/manuscript_figures/figures/manorm_90/annot_", subscript,"_3.tiff", sep=''), 
     # height = 8, width = 12, res=500, units = 'in')
dotplot
dev.off()

head(plot)
df = fread("/Users/trangdo/Google\ Drive/manuscript_figures/figures/manorm_90/enrichment.csv", header=T)
plot = df %>%
  dplyr::mutate(`Fold Enrichment` = as.numeric(`Fold Enrichment`),
                `Gene Ratio` = as.numeric(round(nGenes/`Pathway Genes`, 2))) %>%
  dplyr::arrange(`Fold Enrichment`) %>%
  dplyr::top_n(-25, `Enrichment FDR`) %>%
  dplyr::mutate(Pathway = factor(Pathway, levels=Pathway))  %>%
  ggplot(aes(x = `Fold Enrichment`, y = Pathway, col = -log(`Enrichment FDR`), size = `Gene Ratio`)) + 
  geom_point() +
  geom_segment(aes(xend = 0, yend = Pathway), size = 1) + 
  scale_color_viridis(option = "A", direction = -1, begin = 0.1, end = 0.9,
                      name=expression(paste(-Log[10], FDR, sep=''))) +
  scale_size(name='Gene Ratio') +
  xlim(c(0, 1.9)) +
  theme(aspect.ratio = 2.5,
        axis.title.y = element_blank(),
        axis.title.x = element_text(face='bold', size=12),
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = 'white'),
        legend.position = 'right',
        axis.text.y = element_text(size=12, color = 'black'),
        legend.box="vertical"
        # legend.margin=margin()
        ) 

# plot
tiff(paste("/Users/trangdo/Google\ Drive/manuscript_figures/figures/manorm_90/global_his", subscript,"3.tiff", sep=''), 
     height = 6, width = 9, res=300, units = 'in')
plot
dev.off()
# annot_genes_tissue_cluster = readRDS("annot_genes_tissue_cluster.RDS")
# gene_list_tissue = lapply(annot_genes_tissue_cluster, get_entrez_id)
# names(gene_list_tissue) = histone_type_list
# ck_bp_005_tissue = compareCluster(geneCluster = gene_list_tissue, fun = "enrichGO",
#                                   OrgDb='org.Hs.eg.db', ont = "BP", qvalueCutoff = 0.05,
#                                   pAdjustMethod = "fdr", readable =TRUE)
# saveRDS(ck_bp_005_tissue, "annot_genes_tissue_cluster_bp005.RDS")

#----- Find candidate genes -----
temp = all_pearcor_r$H3K4me3
temp = temp$aorta_CD4positivealphabetaTcell
apply(temp[2:length(temp)], 2, function(x){
  res = !is.na(x) & abs(x) > 0.9
  if (length(res) > 0) {
    print(paste(temp$gene_id[res], x[res]))
  }
})

lapply(all_pearcor_padj_sig, function(x) print(unique(x[x$gene_id == "FGFR2",])))
all_pearcor_padj
grep("FGFR2", unique(unlist(unlist(all_genes_joined_padj))))
"FGFR2" %in% unique(unlist(unlist(all_pearcor_padj_sig)))

sapply(all_genes_joined_padj[3], function(x) lapply(x, function(y) if ("FGFR2" %in% y){print(y)}))
temp = all_pearcor_r[[3]]
temp[temp$gene_id == 'FGFR2', grep(,colnames(temp))]



