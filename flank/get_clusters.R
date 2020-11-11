library(data.table, quietly=TRUE)
library(tidyverse)
library(dplyr, quietly=TRUE)
library(plyr)
library(rtracklayer)
library(ggplot2)
library(reshape2)
library(NMF)
check_subcluster <- function(cluster1, all_cluster_str){
  if (TRUE %in% unique(sapply(all_cluster_str, function(x) all(cluster1 %in% x)))) return(TRUE)
  else return(FALSE)
}
get_adj_mat <- function(all_pairs){
  all_pairs = str_split(all_pairs, ',')[[1]]
  # print(all_pairs)
  all_tissues = Reduce(union, sapply(all_pairs, function(x) str_split(x, '_')))
  # print(all_tissues)
  
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
check_cluster <- function(clusters){
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
      if (!is.null(new_cluster)){
        new_cluster = unique(new_cluster)
        all_cluster_str[[i]] = new_cluster
      }
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

all_genes_clusters = vector("list")
for (k in 1:length(all_genewise_cluster)){
  all_genewise_cluster_H = all_genewise_cluster[[k]]
  all_results = vector("list")
  for (h in 1:length(all_genewise_cluster_H)){
    gene_cluster = all_genewise_cluster_H[h]
    print(names(gene_cluster))
    all_tissues = Reduce(union, sapply(str_split(gene_cluster, ',')[[1]], function(x) str_split(x, '_')))
    adj_mat = get_adj_mat(gene_cluster)
    all_tissues_combi = unlist(Map(combn, list(all_tissues), seq(3, length(all_tissues)), simplify = FALSE), recursive=FALSE)
    all_tissues_combi = all_tissues_combi[order(sapply(all_tissues_combi, length), decreasing=T)]
    gene_cluster_list = check_cluster(all_tissues_combi)
    all_results[[h]] = gene_cluster_list
  }
  all_genes_clusters[[k]] = all_results
}
saveRDS(all_genes_clusters, "all_genes_clusters.RDS")
