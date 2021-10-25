library(data.table, quietly=TRUE)
library(dplyr, quietly=TRUE)
library(boot)
library(stats)
library(parallel)
library(tidyverse)

histone_type_list = c("H3K36me3", "H3K4me3")
check_gene = 'FGFR2'
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
get_all_pairs.his <- function(all_pairs.his, his, check_gene){
  pair.his_list = vector("list", length(all_pairs.his))
  for (i in 1:length(all_pairs.his)){
    # for (i in 1:1){
    print(paste("Pair: ", i, sep=''))
    pair.his = all_pairs.his[[i]]
    pair.his = fread(pair.his)
    print(head(pair.his))
    id = as.data.frame(do.call(rbind, lapply(pair.his$V9, function(x) strsplit(x, split='"', fixed=T)[[1]][c(2, 6)])))
    colnames(id) = c('gene', 'exon')
    pair.his = pair.his %>%
      dplyr::mutate(
        gene = id$gene, exon = id$exon, type = V3,
        p_val = as.numeric(as.character(V14)),
        m_val = as.numeric(as.character(V13))
      ) %>%
      dplyr::select(gene, exon, m_val) %>%
      dplyr::filter(gene == check_gene)
    pair.his_list[[i]] = pair.his
  }
  lapply(pair.his_list, function(x) print(dim(x)))
  return(pair.his_list)
}

get_all_pairs.his_list <- function(histone_type_list, check_gene){
  all_pairs.his_list = vector("list", length(histone_type_list))
  for (j in 1:length(histone_type_list)){
    his = histone_type_list[[j]]
    all_pairs.his = list.files(paste("/home/dhthutrang/ENCODE/chip_seq", his, "flank/fl", sep='/'), pattern = '.txt', full.names = TRUE)
    all_pairs.his_list[[j]] = get_all_pairs.his(all_pairs.his, his, check_gene)
  }
  return(all_pairs.his_list)
}

list_results = get_all_pairs.his_list(histone_type_list, check_gene = 'FGFR2')
saveRDS(list_results, paste(gene, 'manorm.RDS', sep=''))


