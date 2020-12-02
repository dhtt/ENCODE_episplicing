library(data.table, quietly=TRUE)
library(tidyverse)
library(dplyr, quietly=TRUE)
library(plyr)
library(rtracklayer)
library(ggplot2)
library(reshape2)
library(NMF)
library(formattable)
#Annotation
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
library(viridis)
library(ggplot2)
library(ggraph)

library(ggpubr) #plot1 plot2 are in analyze_flank.R 
# setwd("/Users/dhthutrang/Documents/BIOINFO/Episplicing/ENCODE_episplicing/flank")
library(boot)
library(stats)
library(parallel)
library("doMC")
doMC::registerDoMC(cores = 50)
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
# all_res_list.pearcor_padj_sig = readRDS("all_res_list.pearcor_padj_sig.RDS")
# all_genewise_cluster = get_genewise_summary(all_res_list.pearcor_padj_sig)
# saveRDS(all_genewise_cluster, "all_genewise_cluster.RDS")
all_res_list.pearcor_r = readRDS("all_res_list.pearcor_r.RDS")
all_genewise_cluster = get_genewise_summary(all_res_list.pearcor_r)
saveRDS(all_genewise_cluster, "all_genewise_cluster_r.RDS")

all_genewise_cluster = readRDS("all_genewise_cluster_r.RDS")

all_mat_hist = vector("list")
for (k in 1:length(all_genewise_cluster)){
  all_genewise_cluster_H = all_genewise_cluster[[k]]
  all_adj_mat = vector("list")
  for (h in 1:length(all_genewise_cluster_H)){
    gene_cluster = all_genewise_cluster_H[h]
    print(names(gene_cluster))
    all_tissues = Reduce(union, sapply(str_split(gene_cluster, ',')[[1]], function(x) str_split(x, '_')))
    adj_mat = get_adj_mat(gene_cluster, all_tissues)
    all_adj_mat[[h]] = adj_mat
  }
  all_mat_hist[[k]] = all_adj_mat
}
saveRDS(all_mat_hist, "all_mat_hist_r.RDS")
all_mat_hist = readRDS("all_mat_hist_r.RDS")
# print("Length all_mat_hist")
# print(length(all_mat_hist))
# sapply(all_mat_hist, function(x) print(length(x)))

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
saveRDS(all_tissues_hist, "all_tissues_hist_r.RDS")
all_tissues_hist = readRDS("all_tissues_hist_r.RDS")
head(all_tissues_hist[[1]])
print("Length all_tissues_hist")
print(length(all_tissues_hist))
sapply(all_tissues_hist, function(x) print(length(x)))

print("====================================================")
make_cluster_pal <- function(){
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
  for (k in 1:6){
    names(all_genes_clusters[[k]]) = names(all_genewise_cluster[[k]])
  }
}

saveRDS(all_genes_clusters, "all_genes_clusters_pal_named_r.RDS")

all_genes_clusters = readRDS("all_genes_clusters_pal_named_r.RDS")
print(head(all_genes_clusters[[1]]))
n_clusters = sapply(all_genes_clusters[[1]], length)
print(summary(n_clusters))
n_clusters = n_clusters[order(n_clusters, decreasing = TRUE)]
head(n_clusters)
# n_clusters[1:4]

make_upairs <- function(list_tissue){
  all_pairs = c()
  for (tissue1 in list_tissue){
    for (tissue2 in list_tissue){
      if (tissue1 < tissue2) join_pairs = paste(tissue1, tissue2, sep='_')
      else join_pairs=NULL
      all_pairs = c(all_pairs, join_pairs)
    }
  }
  return(all_pairs)
}

official_name = c("Adipose tissue", "Aorta", "CD4-positive alpha beta T cell", 
                  "CD8 positive alpha beta T cell", "Ectodermal cell", "Endodermal cell", 
                  "Esophagus", "H1 cell", "Mesenchymal stem cell", "Mesendoderm", "Mesodermal cell", 
                  "Neuronal stem cell", "Pancreas", "Psoas muscle", "Sigmoid colon", 
                  "Small intestine", "Spleen", "Stomach", "Trophoblast")
epigenomes_annot = read.csv("/home/dhthutrang/ENCODE/utilities/epi_info.csv", 
                            header = TRUE, sep=';', row.names = 1)
rownames(epigenomes_annot)[rownames(epigenomes_annot) == "trophoblastcell"] = "trophoblast"
epigenomes_annot$official_name = official_name

get_genewise_clusters_df <- function(all_genes_clusters, all_genewise_cluster){
  all_genewise_clusters = vector("list")
  all_pairs_id = make_upairs(rownames(epigenomes_annot))
  for (k in 1:6){
    histone_type = histone_type_list[k]
    print(histone_type)
    all_genes_clusters_H = all_genes_clusters[[k]] #for clusters >= 3
    all_genewise_cluster_H = all_genewise_cluster[[k]] #for [pairs]
    
    temp = lapply(all_genes_clusters_H, function(x) unlist(lapply(x, function(y) paste(y, collapse = "_"))) )
    all_clusters = Reduce(union, temp)
    all_clusters = all_clusters[sapply(all_clusters, function(x) length(str_split(x, "_")[[1]])) > 2]
    
    all_gene_for_clusters = sapply(all_clusters, function(x) names(temp[sapply(temp, function(y) x %in% y)]))
    all_gene_for_pairs = sapply( all_pairs_id, function(x) names(all_genewise_cluster_H)[sapply(all_genewise_cluster_H, function(y) length(grep(pattern = x, x=y))) > 0])
    
    all_gene_and_cluster = c(all_gene_for_clusters, all_gene_for_pairs)
    
    all_gene_and_cluster_df = as.data.table(sapply(all_gene_and_cluster, function(x) paste(x, collapse=", ")))
    all_gene_and_cluster_df$clusters = names(all_gene_and_cluster)
    colnames(all_gene_and_cluster_df) = c("genes", "cluster")
    genewise_clusters_df = all_gene_and_cluster_df %>%
      dplyr::select(genes, cluster) %>%
      group_by(genes) %>%
      dplyr::mutate(clusters = paste(cluster, collapse = ', '),
                    n_clusters = n()) %>%
      dplyr::select(-cluster) %>%
      unique() %>%
      ungroup()
    genewise_clusters_df$n_genes = sapply(genewise_clusters_df$genes, function(x) length(str_split(x, ", ")[[1]]))
    genewise_clusters_df$histone_type = histone_type_list[[k]][1]
    
    all_genewise_clusters[[k]] = genewise_clusters_df
  }
  all_genewise_clusters = do.call("rbind", all_genewise_clusters)
  all_genewise_clusters = all_genewise_clusters[all_genewise_clusters$genes != "", ]
  all_genewise_clusters$tissues = lapply(lapply(lapply(all_genewise_clusters$clusters, 
                                                          function(x) str_split(x, ", ")[[1]]), 
                                                   function(y) str_split(y, "_")), 
                                            function(z) unique(unlist(z)))
  all_genewise_clusters$n_tissues = lapply(all_genewise_clusters$tissues, length)
  return(all_genewise_clusters)
}
all_genewise_clusters_df = get_genewise_clusters_df(all_genes_clusters, all_genewise_cluster)
saveRDS(all_genewise_clusters_df, "all_genewise_clusters_df_r.RDS")
all_genewise_clusters_df = readRDS("all_genewise_clusters_df_r.RDS")
all_genewise_clusters_df$n_tissues = unlist(all_genewise_clusters_df$n_tissues)

plot1 = ggplot(data=all_genewise_clusters_df[all_genewise_clusters_df$n_tissues >= 2,], aes(x = n_genes, fill = histone_type)) +
  geom_histogram(binwidth = 1, position = "dodge", alpha = 0.9, col = "black", boundary=0) +
  labs(fill = "Histone mark", title = "Number of gene sets shared between more than 2 tissues") +
  xlab("Size of genes set") + ylab("Occurences") +
  # scale_y_continuous(breaks=seq(0, 1000, 100)) +
  theme_bw() + scale_fill_viridis_d(begin = 0.25, end = 1) +
  theme(aspect.ratio=1, plot.margin	= unit(c(0.2,0,3.15,0), "cm"))
plot1

plot2 = ggplot(data=all_genewise_clusters_df[all_genewise_clusters_df$n_tissues >= 2,], aes(x = n_tissues, fill = histone_type)) +
  geom_histogram(binwidth = 1, position = "dodge", alpha = 0.9, col = "black") +
  labs(fill = "Histone mark", title = "Number of total tissues in clusters sharing the same genes sets") +
  xlab("Number of tissue clusters") + ylab("Occurences") +
  theme_bw() + scale_fill_viridis_d(begin = 0.25, end = 1) +
  # scale_y_continuous(breaks=seq(0, 250, 25)) +
  scale_x_continuous(breaks=seq(0, max(all_genewise_clusters_df$n_tissues), 1)) +
  theme(aspect.ratio=1, plot.margin	= unit(c(0.2,0,3.15,0), "cm")) 
  # geom_vline(data = annot_genes_gene_cluster, xintercept=annot_genes_gene_cluster$min_tissue)
plot2

get_all_tissue_counts <- function(all_genewise_clusters_df){
  tissue_counts = vector("list")
  for (idx in 1:6){
    histone = histone_type_list[[idx]][1]
    tissues_count = all_genewise_clusters_df[all_genewise_clusters_df$histone_type == histone &
                                               all_genewise_clusters_df$n_tissues >= 3, ]
    tissues_count = paste(lapply(tissues_count$clusters, function(x) unlist(x)), collapse = ', ')
    tissues_count = lapply(tissues_count, function(x) gsub('_', ', ', x))
    tissues_count_df = data.frame(tissue = str_split(tissues_count, ', ')[[1]])
    tissues_count_df$histone_type = histone
    tissue_counts[[idx]] = tissues_count_df
  }
  tissue_counts = do.call('rbind', tissue_counts)
  return(tissue_counts)
}

all_tissues_counts = get_all_tissue_counts(all_genewise_clusters_df)
head(all_tissues_counts)

plot4 = ggplot(data=all_tissues_counts, aes(x = tissue , fill = histone_type)) +
  geom_histogram(position = "dodge", alpha = 1, col = "black", stat="count", orientation='x') +
  labs(fill = "Histone mark", title="Occurences of tissue in all tissue clusters") +
  xlab("Tissue") + ylab("Occurences in tissue clusters") +
  theme_bw() + scale_fill_viridis_d(begin = 0, end = 1, option="D")  +
  theme(axis.text.x=element_text(colour="black", size = 10, angle = 45, vjust=0.9, hjust = 1),
        aspect.ratio=1)
plot4

tiff("sup1_r.tiff", width = 18, height = 6, units = "in", res = 200) #save pdf 20*8
figure <- ggarrange(plot1, plot2, plot4, labels = c("A", "B", "C"), ncol = 3)
figure
dev.off()


#------Get genes for annot------
get_entrez_id <-function(gene_list){
  entrez_gene_list = select(org.Hs.eg.db,
                            keys = gene_list,
                            columns = c("ENTREZID", "SYMBOL"),
                            keytype = "SYMBOL")
  return(entrez_gene_list[[2]])
}

#Genes with more than 3 occurences for annotation
annot_genes_gene_cluster = all_genewise_clusters_df %>%
  dplyr::filter(n_tissues > 2) %>%
  dplyr::group_by(histone_type) %>% 
  dplyr::mutate(annot_genes = paste(genes, collapse = ', '),
                annot_genes = gsub('+', ', ', annot_genes, fixed = TRUE)) %>%
  ungroup %>%
  dplyr::select(annot_genes, histone_type) %>%
  unique()

annot_genes_gene_cluster = lapply(annot_genes_gene_cluster$annot_genes, function(x) unique(str_split(x, ', ')[[1]]))
lapply(annot_genes_gene_cluster, length)
paste(annot_genes_gene_cluster[[4]], collapse = ', ')

gene_list_gene = lapply(annot_genes_gene_cluster, get_entrez_id)
sapply(gene_list_gene, function(x) length(x[is.na(x)]))
names(gene_list_gene) = histone_type_list

ck_bp_005_gene = compareCluster(geneCluster = gene_list_gene, fun = "enrichGO",
                                OrgDb='org.Hs.eg.db', ont = "BP", pvalueCutoff = 0.05,
                                pAdjustMethod = "fdr", readable =TRUE)
saveRDS(ck_bp_005_gene, "annot_genes_gene_cluster_bp005_new.RDS")
# ck_bp_0001_gene = compareCluster(geneCluster = gene_list_gene, fun = "enrichGO",
#                                 OrgDb='org.Hs.eg.db', ont = "BP", qvalueCutoff = 0.01,
#                                 pAdjustMethod = "fdr", readable =TRUE)
# saveRDS(ck_bp_0001_gene, "annot_genes_gene_cluster_bp0001_new.RDS")
# ck_bp_0001_gene = readRDS("annot_genes_gene_cluster_bp0001_new.RDS")

plot5 = dotplot(ck_bp_005_gene, showCategory = 25) +
  scale_color_viridis(option = "D") +
  ggtitle(label = "Enriched GO Terms for Tissue-specific Epispliced Genes (Biological Process)") +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x=element_text(colour="black", size = 10, angle = 45, vjust=0.5),
    axis.text.y=element_text(colour="black", size = 10),
    plot.margin = unit(c(20,20,20,20), "pt"))
# plot5

tiff("annot_r.tiff", width = 12, height = 12, units = "in", res = 200) #save pdf 20*8
plot5
dev.off()
