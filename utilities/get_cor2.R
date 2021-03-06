library(data.table, quietly=TRUE)
library(dplyr, quietly=TRUE)
library(boot)
library(stats)
library(parallel)
library("doMC")
#setwd("/Users/dhthutrang/Documents/BIOINFO/Episplicing/ENCODE_episplicing/flank")
doMC::registerDoMC(cores = 17)
#setwd("/home/dhthutrang/ENCODE/flank")

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
all_pairs.exp = list.files("/home/dhthutrang/ENCODE/mRNA_seq/dexseqcount/res", full.names = TRUE)
get_all_pairs.exp <- function(all_pairs.exp){
  colname_exp = c("gene_id", "exon_id", get_colname(all_pairs.exp, "exp"))
  pair.exp_list = vector("list", length(all_pairs.exp))
  for (i in 1:length(all_pairs.exp)){
    print(paste("Pair: ",i, sep=''))
    pair.exp = all_pairs.exp[[i]]
    #fwrite(exp_id, "/home/dhthutrang/ENCODE/utilities/exp_id.txt", col.names = FALSE, row.names = FALSE, quote = FALSE, sep='\t')
    pair.exp_list[[i]] = fread(pair.exp)[,c("stat", "padj")]
  }
  pair.exp_list = lapply(pair.exp_list,
                         function(x) {
                           x = x %>%
                             mutate(exp = if_else(padj <= 0.1 & !is.na(padj), true = stat,
                                                  false = 0.0)) %>%
                             dplyr::select(exp)
                         })
  pair.exp_list = as.data.table(pair.exp_list)
  exp_id = fread("/home/dhthutrang/ENCODE/utilities/exp_id.txt", sep = '\t', quote=FALSE, header = FALSE)
  print(paste("COMPARE LENGTH", dim(exp_id), dim(pair.exp_list), sep=' '))
  pair.exp_list = cbind(exp_id, pair.exp_list)
  pair.exp_list = pair.exp_list[order(pair.exp_list$V1)]
  colnames(pair.exp_list) = colname_exp
  return(pair.exp_list)
}

#all_pairs.exp = get_all_pairs.exp(all_pairs.exp)
#saveRDS(all_pairs.exp, "/home/dhthutrang/ENCODE/flank/all_pairs.exp.RDS")
all_pairs.exp = readRDS("/home/dhthutrang/ENCODE/flank/all_pairs.exp.RDS")
# all_pairs.exp = readRDS("all_pairs.exp.RDS")
print(head(all_pairs.exp))
paste("CONTROL: ", length(unique(all_pairs.exp[all_pairs.exp$H1_mesenchymalstemcell > 0, ]$gene_id)), sep ='')

#===== PREPARE HIS FILE (6 TOTAL) =====
print("===== PREPARE HIS FILE (6 TOTAL) =====")
his_id = fread("flank_id.txt", sep = '\t', quote=FALSE, header = FALSE)
get_all_pairs.his <- function(all_pairs.his){
  pair.his_list = vector("list", length(all_pairs.his))
  for (i in 1:length(all_pairs.his)){
    # for (i in 1:1){
    print(paste("Pair: ", i, sep=''))
    pair.his = all_pairs.his[[i]]
    pair.his = fread(pair.his)
    print(head(pair.his))
    pair.his = pair.his %>%
      mutate(temp_val = abs(as.numeric(as.character(V10))),
	     temp_p = as.numeric(as.character(V11)),
	     m_val = if_else(!is.na(temp_val) & temp_val >= 1 , #& !is.na(temp_p) & temp_p <= 0.1
	                     true = temp_val, false = 0)) %>%
      dplyr::select(m_val)
    pair.his_list[[i]] = pair.his
  }
  pair.his_list = as.data.table(pair.his_list)
  pair.his_list = pair.his_list %>%
    group_by(group = gl(n()/2, 2)) %>%
    summarise_all(max) %>%
    dplyr::select(-group)
  pair.his_list = cbind(his_id, pair.his_list)
  pair.his_list = pair.his_list[order(pair.his_list$V1)]
  return(pair.his_list)
}
get_all_pairs.his_list <- function(histone_type_list){
  all_pairs.his_list = vector("list", length(histone_type_list))
  for (j in 1:length(histone_type_list)){
    # for (j in 1:1){
    his = histone_type_list[[j]]
    print(his)
    all_pairs.his = list.files(paste("/home/dhthutrang/ENCODE/chip_seq", his, "flank", sep='/'), pattern = '.txt', full.names = TRUE)
    print(all_pairs.his)
    colname_his = c("gene_id", "exon_id", get_colname(all_pairs.his, "his")) 
    all_pairs.his.sig = get_all_pairs.his(all_pairs.his)
    colnames(all_pairs.his.sig) = colname_his
    all_pairs.his_list[[j]] = all_pairs.his.sig
  }
  return(all_pairs.his_list)
}
#histone_type_list = list("H3K4me1", "H3K4me3", "H3K9me3", "H3K27me3", "H3K36me3", "H3K27ac")
histone_type_list = list("H3K27ac", "H3K27me3", "H3K36me3", "H3K4me1", "H3K4me3", "H3K9me3")
#all_pairs.his_list = get_all_pairs.his_list(histone_type_list)
#saveRDS(all_pairs.his_list, "/home/dhthutrang/ENCODE/flank/all_pairs.his_list.RDS")
#all_pairs.his_list = readRDS("/home/dhthutrang/ENCODE/flank/all_pairs.his_list.RDS")
all_pairs.his_list = readRDS("all_pairs.his_list.RDS")
head(all_pairs.his_list[[1]], 50)

#===== CORRELATION WITH RANDOMIZATION =====
# ------------ Execute analysis ------------
all_genes = fread("gene_id.txt", header = FALSE)
p_value_calculator <- function(r, nrow){
  P <- r*sqrt(nrow-2)/sqrt(1-r*r)
  P <- 2*pt(-abs(P), nrow-2)
  return(P)
}
pearcor_p <- function(exp, his){
  if (length(unique(exp)) > 1 & length(unique(his)) > 1){
    p_val = p_value_calculator(cor(exp, his, method = "pearson"), nrow = length(exp))
    # p_val = p.adjust(p_val, method = "fdr", n=ncol(all_genes))
    return(p_val)
  }
  else {
    return(NA)
  }
}

pearcor_r <- function(exp, his, n_his){
  df = data.frame(cbind(exp, his))
  n_points = nrow(unique(df))
  if (0 %in% apply(df, 1, unique)) n_points = n_points-1
  #print(n_points)
  if (n_points > 2 & length(unique(exp)) > 1 & length(unique(his)) > n_his){
    r_val = cor(exp, his, method = "pearson")
    # p_val = p.adjust(p_val, method = "fdr", n=ncol(all_genes))
    return(r_val)
  }
  else {
    return(NA)
  }
}

analyze_array <- function(all_pairs.exp, all_pairs.his, n_his){
  all_res_pair = vector("list", ncol(all_pairs.exp) - 2)
  subset_name = colnames(all_pairs.his)
  colnames(all_pairs.exp) = gsub('trophoblastcell', 'trophoblast', colnames(all_pairs.exp))
  all_pairs.exp_subset = all_pairs.exp[, ..subset_name]
  
  #for (i in 1:n_pairs){
  all_res_pair <- foreach( i=1:(length(subset_name)-2), .combine='c', .packages=c('dplyr') ) %dopar% {
    print(paste("Pair: ", i, sep=''))
    exp = all_pairs.exp_subset[[i+2]]
    his = all_pairs.his[[i+2]]
    data_table = as.data.table(cbind(exp, his))
    head(data_table)
    
    res_table = data_table %>%
      group_by(all_pairs.exp$gene_id) %>%
      summarise(res = pearcor_r(exp, his, n_his)) %>%
      dplyr::select(res)
    #all_res_pair[[i]] = res_table
  }
  all_res_pair = as.data.table(all_res_pair)
  all_res_pair = cbind(all_genes, all_res_pair)
  colnames(all_res_pair) = c('gene_id', subset_name[3:length(subset_name)])
  print(head(all_res_pair))
  return(as.data.frame(all_res_pair))
}
analyze_array_list <- function(all_pairs.exp, all_pairs.his_list, n_his){
  all_res_list = vector("list", length(histone_type_list)-1 )
  for (j in 1:length(histone_type_list)){
    print(paste("Histone: ", histone_type_list[[j]], sep = ''))
    all_pairs.his = all_pairs.his_list[[j]]
    all_res_pair = analyze_array(all_pairs.exp, all_pairs.his, n_his)
    all_res_list[[j]] = all_res_pair
  }
  return(all_res_list)
}

print("Pearsons-p correlation")
#all_res_list.pearcor_p = analyze_array_list(all_pairs.exp, all_pairs.his_list, "pearcor_p")
#saveRDS(all_res_list.pearcor_p, "/home/dhthutrang/ENCODE/flank/all_res_list.pearcor_p.RDS")

all_res_list.pearcor_r = analyze_array_list(all_pairs.exp, all_pairs.his_list, 9)
saveRDS(all_res_list.pearcor_r, "/home/dhthutrang/ENCODE/flank/all_res_list.pearcor_r_10.RDS")




