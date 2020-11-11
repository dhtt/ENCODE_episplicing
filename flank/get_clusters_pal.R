library(data.table, quietly=TRUE)
library(tidyverse)
library(dplyr, quietly=TRUE)
library(plyr)
library(boot)
library(stats)
library(parallel)
library("doMC")
#setwd("/Users/dhthutrang/Documents/BIOINFO/Episplicing/ENCODE_episplicing/flank")
doMC::registerDoMC(cores = 17)
histone_type_list = list("H3K27ac", "H3K27me3", "H3K36me3", "H3K4me1", "H3K4me3", "H3K9me3")

check_edge <- function(adj_mat, t1, t2){
  if (t2 %in% adj_mat[[t1]] | t1 %in% adj_mat[[t2]] ) return(TRUE)
  else return(FALSE)
}
check_subcluster <- function(cluster1, all_cluster_str){
  if (TRUE %in% unique(sapply(all_cluster_str, function(x) all(cluster1 %in% x)))) return(TRUE)
  else return(FALSE)
}
get_adj_mat <- function(all_pairs, all_tissues){
  all_pairs = str_split(all_pairs, ',')[[1]]
  
  adj_mat = vector("list", length(all_tissues))
  names(adj_mat) = all_tissues
  for (pair in all_pairs){
    pair = str_split(pair, '_')[[1]]
    t1 = pair[1]
    t2 = pair[2]
    adj_mat[[t1]] = c(adj_mat[[t1]], t2)
    adj_mat[[t2]] = c(adj_mat[[t2]], t1)
  }
  for (i in 1:length(adj_mat)){
    adj_mat_temp = c(adj_mat[[i]], names(adj_mat)[i])
    adj_mat[[i]] = adj_mat_temp[order(adj_mat_temp)]
  }
  adj_mat = adj_mat[order(sapply(adj_mat, length), decreasing=T)]
  return(adj_mat)
}
check_cluster <- function(clusters, adj_mat){
  all_cluster_str = vector("list")
  for (i in 1:length(clusters)){
    cluster = clusters[[i]]
    break_sig = FALSE
    new_cluster = NULL
    
    if (check_subcluster(cluster, all_cluster_str) == FALSE) { #if cluster is not a subset of an existing closed clusters
      for (t1 in cluster){
        for (t2 in cluster){
          if (t2 > t1){
            if (check_edge(adj_mat, t1, t2) == TRUE) new_cluster = c(new_cluster, t1 ,t2)
            else {
              break_sig = TRUE #If there is no edge between any two tissues in cluster -> break cluster iter
              new_cluster = NULL #No new cluster if two tissues are not connected
            } 
            if (break_sig == TRUE) break
          }
          if (break_sig == TRUE) break
        }
        if (break_sig == TRUE) break
      }
    }
    if (!is.null(new_cluster)){
      new_cluster = unique(new_cluster)
      all_cluster_str[[i]] = new_cluster
    }
  }
  all_cluster_str = all_cluster_str[sapply(all_cluster_str, length) > 0]
  return(all_cluster_str)
}

get_genewise_summary <- function(all_genes_joined){
  all_genewise_cluster = vector("list")
  for (idx in (1:6)){
    all_genes_H = Reduce(union, all_genes_joined[[idx]])
    gene_cluster = lapply(all_genes_H, function(y) lapply(all_genes_joined[[idx]], function(x)  y %in% x))
    names(gene_cluster) = all_genes_H
    gene_cluster_col = sapply(gene_cluster, function(x) {
      all_tissues = names(x[x == TRUE])
      return(paste(all_tissues, collapse = ','))})
    
    all_genewise_cluster[[idx]] = gene_cluster_col
  }
  return(all_genewise_cluster)
}
# all_res_list.pearcor_padj_sig = loadRDS("all_res_list.pearcor_padj_sig.RDS")
# all_genewise_cluster = get_genewise_summary(all_res_list.pearcor_padj_sig)
# saveRDS(all_genewise_cluster, "all_genewise_cluster.RDS")
all_genewise_cluster = readRDS("all_genewise_cluster.RDS")
# 
# all_mat_hist = vector("list")
# for (k in 1:length(all_genewise_cluster)){
#   all_genewise_cluster_H = all_genewise_cluster[[k]]
#   all_adj_mat = vector("list")
#   for (h in 1:length(all_genewise_cluster_H)){
#     gene_cluster = all_genewise_cluster_H[h]
#     print(names(gene_cluster))
#     all_tissues = Reduce(union, sapply(str_split(gene_cluster, ',')[[1]], function(x) str_split(x, '_')))
#     adj_mat = get_adj_mat(gene_cluster, all_tissues)
#     all_adj_mat[[h]] = adj_mat
#   }
#   all_mat_hist[[k]] = all_adj_mat
# }
# saveRDS(all_mat_hist, "all_mat_hist.RDS")
all_mat_hist = readRDS("all_mat_hist.RDS")
print("Length all_mat_hist")
print(length(all_mat_hist))
sapply(all_mat_hist, function(x) print(length(x)))

all_tissues_hist = vector("list")
for (k in 1:length(all_genewise_cluster)){
  all_genewise_cluster_H = all_genewise_cluster[[k]]
  all_tissues_H = vector("list")
  for (h in 1:length(all_genewise_cluster_H)){
    gene_cluster = all_genewise_cluster_H[h]
    all_tissues_H[[h]] = Reduce(union, sapply(str_split(gene_cluster, ',')[[1]], function(x) str_split(x, '_')))
  }
  all_tissues_hist[[k]] = all_tissues_H
}
saveRDS(all_tissues_hist, "all_tissues_hist.RDS")
all_tissues_hist = readRDS("all_tissues_hist.RDS")
print("Length all_tissues_hist")
print(length(all_tissues_hist))
sapply(all_tissues_hist, function(x) print(length(x)))

print("====================================================")
all_genes_clusters = vector("list")
# for (k in 1:length(all_genewise_cluster)){
for (k in 1:6){
  print(paste("HISTONE: ", histone_type_list[k], sep=''))
  all_genewise_cluster_H = all_genewise_cluster[[k]]
  all_results <- foreach( h=1:(length(all_genewise_cluster_H)) ) %dopar% {
  # all_results <- foreach( h=1:100 ) %dopar% {
    gene_cluster = all_genewise_cluster_H[h]
    all_tissues = all_tissues_hist[[k]][[h]]
    adj_mat_H = all_mat_hist[[k]][[h]]
    if (length(all_tissues) >= 3) {
      all_tissues_combi = unlist(Map(combn, list(all_tissues), seq(3, length(all_tissues)), simplify = FALSE), recursive=FALSE)
    }
    else {
      all_tissues_combi = unlist(Map(combn, list(all_tissues), seq(2, length(all_tissues)), simplify = FALSE), recursive=FALSE)
    }
    all_tissues_combi = all_tissues_combi[order(sapply(all_tissues_combi, length), decreasing=T)]
    gene_cluster_list = check_cluster(all_tissues_combi, adj_mat_H)
  }
  print("FINALLY FINISHED")
  print(paste("LENGTH RESULT: ", length(all_results), sep=''))
  file_names = paste("all_results_", k, ".RDS", sep='')
  saveRDS(all_results, file_names)
  all_genes_clusters[[k]] = all_results
}

saveRDS(all_genes_clusters, "all_genes_clusters_pal.RDS")
print("EXAMPLE")
print(head(all_genes_clusters[[1]]))

for (k in 1:6){
  names(all_genes_clusters[[k]]) = names(all_genewise_cluster[[k]])
}
saveRDS(all_genes_clusters, "all_genes_clusters_pal_named.RDS")





