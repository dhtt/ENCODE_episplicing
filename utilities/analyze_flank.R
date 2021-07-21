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
setwd("/Users/dhthutrang/Documents/BIOINFO/Episplicing/ENCODE_episplicing/flank")
official_name = c("Adipose tissue", "Aorta", "CD4-positive alpha beta T cell", 
                  "CD8 positive alpha beta T cell", "Ectodermal cell", "Endodermal cell", 
                  "Esophagus", "H1 cell", "Mesenchymal stem cell", "Mesendoderm", "Mesodermal cell", 
                  "Neuronal stem cell", "Pancreas", "Psoas muscle", "Sigmoid colon", 
                  "Small intestine", "Spleen", "Stomach", "Trophoblast")
#============START here============
#------- Read in results ------- 
all_res_list.pearcor_p = readRDS("all_res_list.pearcor_p.RDS")
# get_p_adj <- function(all_res_list.pearcor_p){
#   all_res_list.pearcor_p_adj = vector("list")
#   for (i in 1:length(all_res_list.pearcor_p)){
#     gene_id = all_res_list.pearcor_p[[i]]$gene_id
#     res_list = all_res_list.pearcor_p[[i]]
#     res_list_p = apply(res_list[, 2:ncol(res_list)], 1, p.adjust, method = "fdr")
#     res_list_p = as.data.frame(cbind(gene_id, t(res_list_p)))
#     all_res_list.pearcor_p_adj[[i]] = res_list_p
#   }
#   return(all_res_list.pearcor_p_adj)
#   }
# all_res_list.pearcor_padj = get_p_adj(all_res_list.pearcor_p)
# saveRDS(all_res_list.pearcor_padj, "all_res_list.pearcor_padj.RDS")
# table(all_res_list.pearcor_padj[[1]] <= 0.05)

get_tissue_spec_ref <- function(){
  RPKM = read.csv("57epigenomes.RPKM.pc", row.names=1, sep="")
  head(RPKM)
  # epigenomes = c("E003", "E004", "E005", "E006", "E007", "E011","E012","E013", "E016", "E024", "E053","E054", "E065","E066","E071","E079","E094","E095", "E096", "E098", "E100","E105","E106", "E109","E113") #E022-E027
  epigenomes = c("E065", "E038", "E047", "E012", "E011", "E079","E003","E006", "E004", "E007", "E098","E100", "E106","E109","E113","E094","E005") #E022-E027
  
  RPKM = RPKM[, epigenomes]
  normed_RPKM = as.data.frame(t(apply(RPKM, 1, function(x) (x/max(x)))))
  TSI_RPKM = as.data.frame(apply(normed_RPKM, 1, function(x)((length(normed_RPKM)-sum(as.numeric(x)))/length(normed_RPKM))))
  TSI_RPKM$gene = rownames(RPKM)
  head(TSI_RPKM)
  ENS_list = TSI_RPKM[TSI_RPKM[[1]] >= 0.75, 2]
  
  ENS_gene_list = read.delim("Ensembl_v65.Gencode_v10.ENSG.gene_info.txt", header=FALSE)
  head(ENS_gene_list)
  colnames(ENS_gene_list) = c("ens", "chr", "start", "end", "strand", "feature", "symbol", "name")
  TSI_symbols = ENS_gene_list[ENS_gene_list$ens %in% ENS_list, "symbol"]
  print(length(TSI_symbols))
  return(TSI_symbols)
}
TSI_symbols = get_tissue_spec_ref()
histone_type_list = list("H3K27ac", "H3K27me3", "H3K36me3", "H3K4me1", "H3K4me3", "H3K9me3")

# ----Get significant results -----
get_all_res_list_sig <- function(all_res_list, method, r_sig=0.5, p_sig= 0.05){
  all_res_list_sig = vector("list", length(all_res_list))
  for (i in 1:length(all_res_list)) {
    all_res = all_res_list[[i]]
    
    all_res_sig = vector("list", ncol(all_res)-1)
    for (j in 2:ncol(all_res)){
      all_res.col = as.numeric(all_res[[j]])
      if (method == "pearcor"){
        all_res_sig[[j-1]] = all_res[abs(all_res.col) >= r_sig  & !is.na(all_res.col), 1] 
      }
      else if (method == "pearcor_p" ){
        all_res_sig[[j-1]] = all_res[all_res.col <= p_sig  & !is.na(all_res.col), 1]
      }
    }
    names(all_res_sig) = colnames(all_res)[2:length(colnames(all_res))]
    all_res_list_sig[[i]] = all_res_sig
  }
  return(all_res_list_sig)
}
# all_res_list.pearcor_p_sig = get_all_res_list_sig(all_res_list.pearcor_p, "pearcor_p", p_sig=0.05)
# all_res_list.pearcor_padj_sig = get_all_res_list_sig(all_res_list.pearcor_padj, "pearcor_p", p_sig=0.05)

all_res_list.pearcor_sig = readRDS("new_df/all_res_list.pearcor_sig.RDS")
all_res_list.pearcor_padj_sig = all_res_list.pearcor_sig
all_res_list.pearcor_padj_sig[[1]]$adiposetissue_aorta
paste(all_res_list.pearcor_padj_sig[[6]]$neuronalstemcell_spleen, collapse = ', ')
lapply(all_res_list.pearcor_padj_sig, function(x) paste(x$neuronalstemcell_spleen, collapse = '\',\''))

lapply(all_res_list.pearcor_padj_sig, function(x) paste(x$neuronalstemcell_spleen, collapse = '\',\''))

all_sig_genes = reduce(lapply(all_res_list.pearcor_padj_sig, function(x) reduce(x, union)), union)

# 'SEPTIN9' %in% all_res_list.pearcor_padj_sig[[5]]$neuronalstemcell_spleen

# lapply(all_res_list.pearcor_p_sig, length)
# lapply(all_res_list.pearcor_p_sig, function(x) lapply(x, function(y) length(y)))
# lapply(all_res_list.pearcor_p_sig, function(x) lapply(x, function(y) paste(y, collapse = ', ')))
# temp = sapply(all_res_list.pearcor_padj_sig, function(x) sapply(x, function(y) grep('FGFR2', y)))
# temp1 = unlist(temp[[6]])
# temp1 #his6 ncam1

# ----Tissue spec array-----
get_tissue_spec_array <- function(all_res_list.pearcor_p_sig){
  temp2 = vector("list", length(all_res_list.pearcor_p_sig))
  for (i in 1:length(all_res_list.pearcor_p_sig)){
    temp = names(all_res_list.pearcor_p_sig[[i]])
    temp1 = as.vector(sapply(temp, function(x) strsplit(x, split = '_')[[1]]))
    tissue_name = unique(temp1)
    tissue_index = lapply(tissue_name, function(x) grep(x, temp))
    names(tissue_index) = tissue_name
    temp2[[i]] = tissue_index
  }
  return(temp2)
}
tissue_type_list = get_tissue_spec_array(all_res_list.pearcor_padj_sig)
names(tissue_type_list) = histone_type_list

get_genes_for_tissue <- function(res_list){
  all_genes_joined = vector("list") #List of 6 histone, each has sig genes for 25 tissues
  for (i in 1:length(res_list)){ #for each histone type
    res = res_list[[i]]
    tissue_list = tissue_type_list[[i]]
    all_tissues = lapply(tissue_list, function(x) return(Reduce(union, res[x])))
    names(all_tissues) = names(tissue_list)
    all_genes_joined[[i]] = all_tissues
  }
  names(all_genes_joined) = histone_type_list
  return(all_genes_joined)
}
# all_genes_joined = get_genes_for_tissue(all_res_list.pearcor_p_sig)
all_genes_joined_padj = get_genes_for_tissue(all_res_list.pearcor_padj_sig)
# length(all_genes_joined)
# lapply(all_genes_joined, length)
# lapply(all_genes_joined[[5]], length)
# lapply(all_genes_joined, function(x) lapply(x, function(y) length(y)))
# length(all_res_list.pearcor_padj_sig[[1]]$adiposetissue_aorta)

get_all_len_before <- function(all_genes_joined){
  all_len_before = lapply(all_genes_joined, function(x) lapply(x, function(y) length(y)))
  all_len_before = lapply(all_len_before, function(x) as.data.frame(cbind(names(x), do.call(rbind, x))))
  all_len_before = all_len_before %>% purrr::reduce(full_join, by = "V1")
  colnames(all_len_before) = c('Tissue', histone_type_list)
  return(all_len_before)
}
# all_len_before = get_all_len_before(all_genes_joined)
all_len_before_padj = get_all_len_before(all_genes_joined_padj)
all_len_before_padj$Tissue = official_name
# saveRDS(all_len_before_padj, "all_len_before_padj.RDS")
# fwrite(all_len_before_padj, "all_len_before_padj.csv", quote=FALSE, row.names = FALSE, col.names = TRUE, sep='\t')

get_genes_for_tissue_filtered <- function(all_genes_joined){
  for (i in 1:length(all_genes_joined)){
    print(paste("New histone type", histone_type_list[i], sep=' '))
    for (j in 1: length(all_genes_joined[[i]])){
      print(paste(i, j, sep=' '))
      gene_set = all_genes_joined[[i]][[j]]
      gene_set = gene_set[gene_set %in% TSI_symbols]
      all_genes_joined[[i]][[j]] = gene_set
      print(length(gene_set))
    }
    # names(all_genes_joined) = names(all_genes_joined)
  } #List of 6 histone, each has FILTERED sig genes for 25 tissues
  return(all_genes_joined)
}
# all_genes_joined_filtered = get_genes_for_tissue_filtered(all_genes_joined)
all_genes_joined_filtered_padj = get_genes_for_tissue_filtered(all_genes_joined_padj)
# length(all_genes_joined) #gene names
# lapply(all_genes_joined, length)
# lapply(all_genes_joined[[6]], length) 

get_all_len_after <- function(all_genes_joined){
 
  all_len_after = lapply(all_genes_joined, function(x) lapply(x, function(y) length(y)))
  all_len_after = lapply(all_len_after, function(x) as.data.frame(cbind(names(x), do.call(rbind, x))))
  all_len_after = all_len_after %>% purrr::reduce(full_join, by = "V1")
  colnames(all_len_after) = c('Tissue', histone_type_list)
  return(all_len_after)
}
# all_len_after = get_all_len_after(all_genes_joined_filtered)
all_len_after_padj = get_all_len_after(all_genes_joined_filtered_padj)

# ----- Summ after overlap -----
temp = all_len_after_padj[2:7]
temp <- mutate_all(temp, function(x) as.numeric(as.character(x)))
all_len_after_padj[2:7] = temp
all_len_after_padj$`Total epispliced genes (with overlaps)` = rowSums(all_len_after_padj[2:7], na.rm = TRUE)

length(all_genes_joined_filtered_padj)
names(all_genes_joined_filtered_padj[[1]])

epigenomes_annot = read.csv("/Users/dhthutrang/Documents/BIOINFO/Episplicing/ENCODE_episplicing/utilities/epi_info.csv", 
                            header = TRUE, sep=';', row.names = 1)
rownames(epigenomes_annot)[rownames(epigenomes_annot) == "trophoblastcell"] = "trophoblast"
epigenomes_annot$official_name = official_name
epigenomes_annot$official_name

epigenomes_names = rownames(epigenomes_annot)
epigenomes_names[19] = "trophoblast"
all_res_list_sig_joined = vector("list")
for (i in 1:length(epigenomes_names)){
  e_id = epigenomes_names[[i]]
  all_res_sig = vector("list")
  for (j in 1:6){
    e_idx = match(e_id, names(all_genes_joined_filtered_padj[[j]]))
    print( paste(j, i, e_idx, e_id, sep = ", "))
    all_res_sig[[j]] = all_genes_joined_filtered_padj[[j]][[e_idx]]
  }
  temp = Reduce(union, all_res_sig)
  all_res_list_sig_joined[[i]] = temp
}
all_len_after_padj$`Total epispliced genes (without overlaps)` = sapply(all_res_list_sig_joined, length)
all_len_after_padj$Tissue = official_name
saveRDS(all_len_after_padj, "all_len_after_padj.RDS")
fwrite(all_len_after_padj, "all_len_after_padj.csv", quote=FALSE, row.names = FALSE, col.names = TRUE, sep='\t')

# all_len_after_padj = all_len_after_padj[2:nrow(all_len_after_padj),]
# rownames(all_len_after_padj) = NULL
# 
# all_len_before_padj = all_len_before_padj[2:nrow(all_len_before_padj),]
# rownames(all_len_before_padj) = NULL

# -----FORMATTABLE TABLE-----
head(all_len_before_padj)
tissue_formatter <- formatter("span", 
                              style = x ~ style(
                                width = suffix(x, "px"),
                                font.weight = "bold", 
                                color = ifelse(x == "Total", "black", "gray")))
all_len_after_padj$`Total epispliced genes (with overlaps)`[1] = 1369
all_len_after_padj$`Total epispliced genes (without overlaps)`[1] = 1369
formattable(all_len_before_padj,
            # align =c("l", "l", "c","c","c","c"), 
            list(
              `H3K4me1` = color_tile("white", "wheat"),
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






#============END here============

# ----- Heatmap of tissue-spec genes -----
length(intersect(all_genes_joined_filtered[[1]][[1]], all_genes_joined_filtered[[1]][[2]]))
Reduce(intersect, all_genes_joined_filtered[[5]])

all_sims = vector("list")
for (i in 1:length(all_genes_joined_filtered_padj)){
  print(histone_type_list[i])
  temp = all_genes_joined_filtered_padj[[i]]
  distances = sapply(temp, function(x) sapply(temp, function(y) length(intersect(x, y))))
  distances_union = sapply(temp, function(x) sapply(temp, function(y) length(union(x, y))))
  if (i==1) print(distances)
  diag(distances) = NA
  # distances = distances/max(distances, na.rm = TRUE) #initially
  distances = distances/distances_union #Jaccard
  
  # print(names(temp))
  rownames(distances) = sapply(names(temp), function(x) epigenomes_annot$official_name[match(x, rownames(epigenomes_annot))])
  colnames(distances) = rownames(distances)
  
  # print(rownames(distances))
  all_sims[[i]] = distances
}

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
breaks = seq(0, 0.34, 0.0034)
color = '-RdYlBu2:101'
tiff("heatmap_jac_check.tiff", width = 15, height = 15, units="in", res=200)
par(mfrow = c(2,3), mar = c(1,1,1,1)*10)
for (i in 1:6){
  if (i %in% c(1,2,5,6)){
    aheatmap(all_sims[[i]], Rowv = FALSE, Colv = TRUE, scale="none", 
             main=paste(histone_type_list[[i]]), 
             annCol = epigenomes_annot[names(rownames(all_sims[[i]])), c("Potency", "Type", "Origin", "Life.stage")],
             # annRow = epigenomes_annot[names(rownames(all_sims[[i]])), c("Potency", "Type", "Origin", "Life.stage")],
             annColors = epigenomes_colors[[1]], breaks = breaks, color = color, 
             legend=FALSE, annLegend=FALSE, fontsize=14, cexRow=1, cexCol = 1)
  } 
  else {
    aheatmap(all_sims[[i]], Rowv = FALSE, Colv = TRUE, scale="none", 
             main=paste(histone_type_list[[i]]), 
             annCol = epigenomes_annot[names(rownames(all_sims[[i]])), c("Potency", "Type", "Origin", "Life.stage")],
             # annRow = epigenomes_annot[names(rownames(all_sims[[i]])), c("Potency", "Type", "Origin", "Life.stage")],
             annColors = epigenomes_colors[[2]], breaks = breaks, color = color, 
             legend=FALSE, annLegend=FALSE, fontsize=14, cexRow=1, cexCol = 1)
  }
}
dev.off()

# separate images 
folder = "/Users/dhthutrang/Documents/BIOINFO/Episplicing/ENCODE_episplicing/fig_tab"
tiff(paste(paste(folder, histone_type_list[[1]], sep='/'), "tiff", sep='.'), width = 15, height = 12, units="in", res=100)
heatmap1 = aheatmap(all_distances[[1]], Rowv = TRUE, Colv = TRUE, scale="none", 
         main=paste(histone_type_list[[1]]), 
         # annCol = epigenomes_annot[rownames(all_distances[[1]]), ], 
         # annRow = epigenomes_annot[rownames(all_distances[[1]]), ],
         annCol = epigenomes_annot[epigenomes_annot$official_name == rownames(all_distances[[1]]), c("Potency", "Type", "Origin", "Life.stage")], 
         annRow = epigenomes_annot[epigenomes_annot$official_name == rownames(all_distances[[1]]), c("Potency", "Type", "Origin", "Life.stage")],
         annColors = epigenomes_colors[[1]], breaks = breaks, color = color, 
         legend=FALSE, annLegend=FALSE, fontsize=14, cexRow=1, cexCol = 1)
dev.off()
tiff(paste(paste(folder, histone_type_list[[2]], sep='/'), "tiff", sep='.'), width = 12, height = 12, units="in", res=100)
aheatmap(all_distances[[2]], Rowv = TRUE, Colv = TRUE, scale="none", 
         main = paste(histone_type_list[[2]]), annCol = epigenomes_annot, annRow = epigenomes_annot,
         annColors = epigenomes_colors, breaks = breaks, color = color, 
         legend = FALSE, annLegend = FALSE, labRow = NA, cellwidth = 25, cellheight = 15, fontsize = 20)

dev.off()
tiff(paste(paste(folder, histone_type_list[[3]], sep='/'), "tiff", sep='.'), width = 12, height = 12, units="in", res=100)
aheatmap(all_distances[[3]], Rowv = TRUE, Colv = TRUE, scale="none", 
         main = paste(histone_type_list[[3]]), annCol = epigenomes_annot, annRow = epigenomes_annot,
         annColors = epigenomes_colors, breaks = breaks, color = color, 
         legend = FALSE, annLegend = FALSE, labRow = NA, cellwidth = 25, cellheight = 15, fontsize = 20)

dev.off()
tiff(paste(paste(folder, histone_type_list[[4]], sep='/'), "tiff", sep='.'), width = 12, height = 12, units="in", res=100)
aheatmap(all_distances[[4]], Rowv = TRUE, Colv = TRUE, scale="none", 
         main = paste(histone_type_list[[4]]), annCol = epigenomes_annot, annRow = epigenomes_annot,
         annColors = epigenomes_colors, breaks = breaks, color = color, 
         legend = FALSE, annLegend = FALSE, labRow = NA, cellwidth = 25, cellheight = 15, fontsize = 20)

dev.off()
tiff(paste(paste(folder, histone_type_list[[5]], sep='/'), "tiff", sep='.'), width = 12, height = 12, units="in", res=100)
aheatmap(all_distances[[5]], Rowv = TRUE, Colv = TRUE, scale="none", 
         main = paste(histone_type_list[[5]]), annCol = epigenomes_annot, annRow = epigenomes_annot,
         annColors = epigenomes_colors, breaks = breaks, color = color, 
         legend = FALSE, annLegend = FALSE, labRow = NA, cellwidth = 25, cellheight = 15, fontsize = 20)

dev.off()
tiff(paste(paste(folder, histone_type_list[[6]], sep='/'), "tiff", sep='.'), width = 12, height = 12, units="in", res=100)
aheatmap(all_distances[[6]], Rowv = TRUE, Colv = TRUE, scale="none", 
         main=paste(histone_type_list[[6]]), 
         annCol = epigenomes_annot[-c(10,11,12),], annRow = epigenomes_annot[-c(10,11,12),],
         annColors = list(epigenomes_colors[[1]], epigenomes_colors[[2]][-c(7)]),
         breaks = breaks, color = color, 
         legend = FALSE, annLegend = FALSE, labRow = NA, cellwidth = 25, cellheight = 15, fontsize = 20)

dev.off()
tiff(paste(paste(folder, histone_type_list[[7]], sep='/'), "tiff", sep='.'), width = 12, height = 12, units="in", res=100)
aheatmap(all_distances[[7]], Rowv = TRUE, Colv = TRUE, scale="none", 
         main = paste(histone_type_list[[7]]), annCol = epigenomes_annot, annRow = epigenomes_annot,
         annColors = epigenomes_colors, breaks = breaks, color = color, 
         legend = FALSE, annLegend = FALSE, labRow = NA, cellwidth = 25, cellheight = 15, fontsize = 20)

dev.off()

# ----- Annotation of common genes from HEATMAP -----
epigenomes_names[c(2,8,7,3)] #H1 1
epigenomes_names[c(9, 11, 12)] #H1 2
epigenomes_names[c(11, 12, 9, 10)] #H2
epigenomes_names[c(14,13,18)] #H2 2
epigenomes_names[c(16, 20, 25, 19)] #H3 1
epigenomes_names[c(11, 12, 9, 10)] #H3 2
epigenomes_names[c(11, 12, 9, 10, 1)] #H4 1
epigenomes_names[c(3, 4, 6, 2)] #H4 2
epigenomes_names[c(16,23)]
epigenomes_names[c(8,3,7,6)] #H5
epigenomes_names[c(14,25,18)] #H6 1
epigenomes_names[c(24,23,16,19)] #H6 1
epigenomes_names[c(3, 4, 6, 2)] #H6 1
epigenomes_names[c(22, 21, 13,25,24, 20,14, 18, 15)] #H6 1

# cluster = list(list(1, c(2, 8, 7, 3), c(23, 22, 25, 17, 13, 15, 16, 24), c(6, 5, 4), c(21, 19, 20, 18)),
#                list(2, c(18, 23,17,25), c(5, 2, 4), c(24, 20, 19, 16, 22, 21)),
#                list(3, c(23, 13, 17, 15, 22, 21, 5, 24), c(16, 20, 25, 19), c(4, 2, 6), c(1, 12, 10, 11, 9)),
#                list(4, c(1, 12, 10, 11), c(16, 19, 17, 23), c(3, 4, 6, 2), c(22, 21, 13, 25, 24, 20,14, 18, 15)),
#                list(5, c(22, 21, 23), c(8,3,7,6), c(25, 15, 19, 17), c(24, 13, 18)),
#                list(6, c(14, 18, 25), c(5, 8, 7), c(17, 13, 24, 23, 16, 19), c(20, 22, 21)), 
#                list(7, c(20, 23, 21, 18), c(2, 3, 5, 4), c(22, 15, 17, 16), c(14, 24, 25, 19, 13)))
epigenomes_annot
# cluster = list(list(1, c("trophoblastcell", "H1", "neuronalstemcell", "ectodermalcell")),
#                list(2, c("stomach", "sigmoidcolon", "smallintestine", "esophagus", "pancreas"), c("mesenchymalstemcell", "CD4positivealphabetaTcell", "spleen", "CD8positivealphabetaTcell")),
#                list(3, c("stomach", "smallintestine", "sigmoidcolon", "esophagus"), c("CD4positivealphabetaTcell", "CD8positivealphabetaTcell", "mesenchymalstemcell")),
#                list(4, c("pancreas", "stomach", "sigmoidcolon", "smallintestine", "esophagus"), c("mesenchymalstemcell", "CD4positivealphabetaTcell", "spleen"), c("H1", "trophoblastcell", "neuronalstemcell", "mesendoderm")),
#                list(5, c("stomach", "smallintestine", "sigmoidcolon"), c("H1", "neuronalstemcell", "mesendoderm"), c("mesodermalcell", "mesenchymalstemcell", "spleen", "CD4positivealphabetaTcell")),
#                list(6, c("aorta", "psoasmuscle", "mesodermalcell"), c("stomach", "sigmoidcolon", "esophagus", "pancreas"), c("CD4positivealphabetaTcell", "spleen", "mesenchymalstemcell")))

cluster = list(list(1, c("adiposetissue", "smallintestine", "aorta", "stomach", "psoasmuscle", "sigmoidcolon", "esophagus"), c("trophoblastcell", "neuronalstemcell", "mesendoderm", "ectodermalcell")),
               list(2, c("psoasmuscle", "stomach", "sigmoidcolon", "aorta", "smallintestine", "esophagus", "pancreas"), c("neuronalstemcell", "ectodermalcell", "mesenchymalstemcell", "endodermalcell")),
               list(3, c("stomach", "aorta", "smallintestine", "sigmoidcolon", "psoasmuscle", "esophagus"), c("psoasmuscle", "neuronalstemcell", "mesodermalcell", "mesendoderm")),
               list(4, c("pancreas", "stomach", "aorta", "sigmoidcolon", "smallintestine", "esophagus", "psoasmuscle"), c("trophoblastcell", "endodermalcell", "neuronalstemcell", "mesendoderm")),
               list(5, c("stomach", "smallintestine", "aorta", "sigmoidcolon", "psoasmuscle"), c("neuronalstemcell", "mesendoderm", "ectodermalcell", "endodermalcell")),
               list(6, c("smallintestine", "aorta", "psoasmuscle", "stomach", "sigmoidcolon", "esophagus", "pancreas"), c("mesenchymalstemcell", "trophoblastcell", "ectodermalcell", "endodermalcell")))

get_heatmap_subset <- function(joined_filtered, reduce_method = "union"){
  genes_all_clusters = vector("list", length(cluster))
  for (j in 1:length(cluster)){
    histone_type = cluster[[j]][[1]]
    print(histone_type_list[histone_type])
    
    gene_each_cluster = vector("list", length(cluster[[j]]) - 1)
    for (i in 2:length(cluster[[j]])){
      group_name = cluster[[j]][[i]]
      # group_name = paste(epigenomes_names[group], collapse = ', ')
      
      if (reduce_method == "union"){
        subset = Reduce(union, joined_filtered[[histone_type]][group_name])
      }
      else if (reduce_method == "intersect"){
        subset = Reduce(intersect, joined_filtered[[histone_type]][group_name])
      }
      print(group_name)
      print(length(subset))
      gene_each_cluster[[i-1]] = subset
    }
    names(gene_each_cluster) = paste("cluster", seq(1:length(gene_each_cluster)), sep = '_')
    genes_all_clusters[[j]] = gene_each_cluster
  }
  return(genes_all_clusters)
}
genes_all_clusters = get_heatmap_subset(all_genes_joined_filtered, "union")
genes_all_clusters_padj = get_heatmap_subset(all_genes_joined_filtered_padj, "union")

head(genes_all_clusters[[1]])
lapply(genes_all_clusters_padj[[4]], function(x) paste(x, collapse = ', '))
lapply(genes_all_clusters, function(y) lapply(y, function(x) paste(x, collapse = ', ')))


lapply(cluster, function(x) paste(length(Reduce(union, all_genes_joined_filtered[[1]][x]))))
lapply(cluster, function(x) paste(Reduce(intersect, all_genes_joined[[1]][x]), collapse = ', '))

temp = Reduce(intersect, all_genes_joined_filtered[[5]][c(8,3,7,6)])
paste(temp, collapse = ", ")

Reduce(intersect, all_res_list.pearcor_p_sig_joinedtissue[[5]])

#in get_cluster_pal.R

all_genes_H1 = Reduce(union, all_genes_joined_filtered_padj[[1]])

#Get histogram of number of genes share by n tissues/cells
get_genewise_summary <- function(all_genes_joined){
  all_genewise_cluster = vector("list")
  for (idx in (1:6)){
    all_genes_H1 = Reduce(union, all_genes_joined[[idx]])
    gene_cluster = lapply(all_genes_H1, function(y) lapply(all_genes_joined[[idx]], function(x)  y %in% x))
    names(gene_cluster) = all_genes_H1
    gene_cluster = as.data.table(unlist(lapply(gene_cluster, function(x) { all_tissues = names(x[x == TRUE])
                    return(paste(all_tissues, collapse = ', ' ))})))
    gene_cluster$gene = all_genes_H1
    gene_cluster = gene_cluster %>%
      dplyr::select(V1, gene) %>%
      group_by(V1) %>%
      dplyr::mutate(cluster = V1,
                    genes = paste(gene, collapse = ','),
                    n_genes = n()) %>%
      dplyr::select(-gene) %>%
      unique() %>%
      ungroup() %>%
      dplyr::select(-V1)
    gene_cluster$n_tissues = sapply(gene_cluster$cluster, function(x) length(str_split(x, ', ')[[1]]))
    gene_cluster$histone_type = histone_type_list[[idx]][1]
    all_genewise_cluster[[idx]] = gene_cluster
  }
  all_genewise_cluster = do.call("rbind", all_genewise_cluster)
  return(all_genewise_cluster)
}
all_genewise_cluster = get_genewise_summary(all_genes_joined_filtered_padj)
all_genewise_cluster[all_genewise_cluster$n_genes == 19, ]
all_genewise_cluster[all_genewise_cluster$n_tissues == 19, "genes"]
# all_genewise_cluster_unfiltered = get_genewise_summary(all_genes_joined_padj) #only for STIM1/FGFR2 cuz not tissue-spec 
# all_genewise_cluster_unfiltered[grep("FGFR2", all_genewise_cluster_unfiltered$genes), ]
# min(all_genewise_cluster$n_tissues)

all_genewise_cluster[all_genewise_cluster$n_tissues == 2, ]

plot1 = ggplot(data=all_genewise_cluster, aes(x = n_genes, fill = histone_type)) + 
  geom_histogram(binwidth = 1, position = "dodge", alpha = 0.9, col = "black", boundary=0) +  
  labs(fill = "Histone mark") + xlab("Gene cluster size") + ylab("Occurences") +
  geom_vline(xintercept = 3, linetype="dashed") +
  theme_bw() + scale_fill_viridis_d(begin = 0.25, end = 1) +
  expand_limits(x = 1)
plot2 = ggplot(data=all_genewise_cluster, aes(x = n_tissues, fill = histone_type)) + 
  geom_histogram(binwidth = 1, position = "dodge", alpha = 0.9, col = "black", boundary=0) +  
  labs(fill = "Histone mark") + xlab("Epigenome cluster size") + ylab("Occurences") +
  geom_vline(xintercept = 9, linetype="dashed") +
  theme_bw() + scale_fill_viridis_d(begin = 0.25, end = 1) +
  scale_x_continuous(breaks=seq(0, 20, 1)) +
  expand_limits(x = 1)
plot1
plot2

#------Get genes for annot------
# Tissues cluster with more than 3 occurences
annot_genes_tissue_cluster = temp %>%
  dplyr::filter(n_tissues >= 9) %>%
  dplyr::group_by(histone_type) %>%
  dplyr::mutate(annot_genes = paste(genes, collapse = ',')) %>%
  ungroup %>%
  dplyr::select(annot_genes, histone_type) %>%
  unique() 
head(annot_genes_tissue_cluster)
annot_genes_tissue_cluster$annot_genes[1]
annot_genes_tissue_cluster = lapply(annot_genes_tissue_cluster$annot_genes, function(x) str_split(x, ',')[[1]])
saveRDS(annot_genes_tissue_cluster, "annot_genes_tissue_cluster.RDS")
annot_genes_tissue_cluster = readRDS("annot_genes_tissue_cluster.RDS")
annot_genes_tissue_cluster[[6]]
paste(annot_genes_tissue_cluster[[1]], collapse = as.character(", ") )

#Genes with more than 3 occurences for annotation
annot_genes_gene_cluster = temp %>%
  dplyr::filter(n_genes >= 3) %>%
  dplyr::group_by(histone_type) %>%
  dplyr::mutate(annot_genes = paste(genes, collapse = ',')) %>%
  ungroup %>%
  dplyr::select(annot_genes, histone_type) %>%
  unique() 
head(annot_genes_gene_cluster)
annot_genes_gene_cluster$annot_genes[6]
annot_genes_gene_cluster = lapply(annot_genes_gene_cluster$annot_genes, function(x) str_split(x, ',')[[1]])
saveRDS(annot_genes_gene_cluster, "annot_genes_gene_cluster.RDS")
annot_genes_gene_cluster = readRDS("annot_genes_gene_cluster.RDS")


temp4 = temp3[temp3$n_genes >= 2, "genes"]
temp5 = lapply(temp4, function(x) paste(x, collapse = ","))
#------Get annot for genes------
get_entrez_id <-function(gene_list){
  entrez_gene_list = select(org.Hs.eg.db,
                            keys = gene_list,
                            columns = c("ENTREZID", "SYMBOL"),
                            keytype = "SYMBOL")
  return(entrez_gene_list[[2]])
}
prepare_ck <- function(gene_list, option, q_val = 0.05){
  print( names(gene_list))
  gene_list = lapply(gene_list, get_entrez_id)
  # names(gene_list) = c("ESC", "Mes", "Tro", "Gas", "Ven", "Mus")
  # names(gene_list) = gene_list_names
  if (option == "BP"){
    print("Preparing BP cluster")
    ck_bp = compareCluster(geneCluster = gene_list, fun = "enrichGO",
                           OrgDb='org.Hs.eg.db', ont = "BP", qvalueCutoff = q_val,
                           pAdjustMethod = "fdr", readable =TRUE)
    return(ck_bp)
  }
  else if (option == "MF"){
    print("Preparing MF cluster")
    ck_mf <- compareCluster(geneCluster = gene_list, fun = "enrichGO",
                            OrgDb='org.Hs.eg.db', ont = "MF", qvalueCutoff = q_val,
                            pAdjustMethod = "fdr", readable =TRUE)
    return(ck_mf)
  }
  else if (option == "CC"){
    print("Preparing CC cluster")
    ck_cc <- compareCluster(geneCluster = gene_list, fun = "enrichGO",
                            OrgDb='org.Hs.eg.db', ont = "CC", qvalueCutoff = q_val,
                            pAdjustMethod = "fdr", readable =TRUE)
    return(ck_cc)
  }
}

annot_genes_gene_cluster = readRDS("annot_genes_gene_cluster.RDS")
gene_list_gene = lapply(annot_genes_gene_cluster, get_entrez_id)
names(gene_list_gene) = histone_type_list
ck_bp_005_gene = compareCluster(geneCluster = gene_list_gene, fun = "enrichGO",
                                OrgDb='org.Hs.eg.db', ont = "BP", qvalueCutoff = 0.05,
                                pAdjustMethod = "fdr", readable =TRUE)
saveRDS(ck_bp_005_gene, "annot_genes_gene_cluster_bp005.RDS")


annot_genes_tissue_cluster = readRDS("annot_genes_tissue_cluster.RDS")
gene_list_tissue = lapply(annot_genes_tissue_cluster, get_entrez_id)
names(gene_list_tissue) = histone_type_list
ck_bp_005_tissue = compareCluster(geneCluster = gene_list_tissue, fun = "enrichGO",
                           OrgDb='org.Hs.eg.db', ont = "BP", qvalueCutoff = 0.05,
                           pAdjustMethod = "fdr", readable =TRUE)
saveRDS(ck_bp_005_tissue, "annot_genes_tissue_cluster_bp005.RDS")

plot3 = dotplot(ck_bp_005_gene, showCategory = 8) +
  scale_color_viridis(option = "D") +
  # ggtitle("Enriched GO Terms for Tissue-specific Epispliced Genes (Biological Process)") +
  theme(
    plot.title = element_text(size = 11, hjust = 0.5, face = "bold"),
    axis.text.x=element_text(colour="black", size = 10, angle = 45, vjust=0.5),
    axis.text.y=element_text(colour="black", size = 10),
    plot.margin = unit(c(20,20,20,20), "pt"))
plot3

plot4 = dotplot(ck_bp_005_tissue, showCategory = 10) +
  scale_color_viridis(option = "D") +
  # ggtitle("Enriched GO Terms for Tissue-specific Epispliced Genes (Biological Process)") +
  theme(
    plot.title = element_text(size = 11, hjust = 0.5, face = "bold"),
    axis.text.x=element_text(colour="black", size = 10, angle = 45, vjust=0.5),
    axis.text.y=element_text(colour="black", size = 10),
    plot.margin = unit(c(20,20,20,20), "pt"))
plot4

get_all_tissue_counts <- function(){
  tissue_counts = vector("list")
  for (idx in 1:6){
    histone = histone_type_list[[idx]][1]
    tissues_count = all_genewise_cluster[all_genewise_cluster$histone_type == histone &
                                           all_genewise_cluster$n_tissues >= 3, ]
    tissues_count = paste(tissues_count$cluster, collapse = ', ')
    tissues_count_df = data.frame(tissue = str_split(tissues_count, ', ')[[1]])
    tissues_count_df$histone_type = histone
    tissue_counts[[idx]] = tissues_count_df
  }
  tissue_counts = do.call('rbind', tissue_counts)
  return(tissue_counts)
}
all_tissues_counts = get_all_tissue_counts()
head(all_tissues_counts)

plot5 = ggplot(data=all_tissues_counts, aes(x =  tissue, fill = histone_type)) +
  geom_histogram(position = "stack", alpha = 1, col = "black", stat="count", orientation='x') +
  labs(fill = "Histone mark") + xlab("Tissue") + ylab("Occurences in tissue clusters") +
  theme_bw() + scale_fill_viridis_d(begin = 0, end = 1, option="D")  +
  theme(axis.text.x=element_text(colour="black", size = 10, angle = 45, vjust=0.9, hjust = 1))
plot5

tiff("annot.tiff", width = 20, height = 18, units = "in", res = 200) #save pdf 20*8
figure <- ggarrange(ggarrange(plot1, plot2, plot5, labels = c("A", "B", "C"), ncol = 3),
                    ggarrange(plot3, plot4, labels = c("D", "E"), ncol = 2),
                    nrow = 2, heights = c(2,5))
figure
dev.off()


#------Silhoulette value-------


all_distances = vector("list")
for (i in 1:length(all_sims)){
  all_distances[[i]] = 1-all_sims[[i]]
}

clusters = list( )
temp = list(epigenomes_annot$official_name[c(16,2,1,15,13,14,18, 19,7,17,3,4,5,6,8,9,10,11,12)],
            epigenomes_annot$official_name[c(2,15,19,18,14,16,17,13,7,4,3,10,5,9,8,6,12)],
            epigenomes_annot$official_name[c(2,14,15,7,3,4,19,12,6,8,13,18,16,10,9,17,11)],
            epigenomes_annot$official_name[c(2,7,16,18,15,14,12,10,8,9,13,3,17,19,11,6)],
            epigenomes_annot$official_name[c(17,9,14,7,18,15,2,16,13,4,3,8,10,11,5,6,12)],
            epigenomes_annot$official_name[c(2,16,15,7,5,6,11,8,3,4,9,18,13,19,14,17,12)])


H1 = list(list(temp[[1]][1:7], temp[[1]][8], temp[[1]][9:12], temp[[1]][13:19]),
          list(temp[[1]][1], temp[[1]][2:3], temp[[1]][4], temp[[1]][5], temp[[1]][6], temp[[1]][7], temp[[1]][8], temp[[1]][9], temp[[1]][10], temp[[1]][11:12], temp[[1]][13:14],temp[[1]][15],temp[[1]][16],temp[[1]][17],temp[[1]][18],temp[[1]][19]),
          list(temp[[1]][1:7], temp[[1]][8], temp[[1]][9:10], temp[[1]][11:12], temp[[1]][13], temp[[1]][14], temp[[1]][15:19]),
          list(temp[[1]][1:7], temp[[1]][8], temp[[1]][9:12], temp[[1]][13], temp[[1]][14], temp[[1]][15:16]))
H1_s = lapply(H1, function(x) get_s_tissue(x, 1))

H2 = list(list(temp[[2]][1:2], temp[[2]][3], temp[[2]][4:11], temp[[2]][12:17]),
          list(temp[[2]][1], temp[[2]][2], temp[[2]][3], temp[[2]][4], temp[[2]][5], temp[[2]][6], temp[[2]][7], temp[[2]][8:9], temp[[2]][10:11], temp[[2]][12], temp[[2]][13], temp[[2]][14], temp[[2]][15], temp[[2]][16], temp[[2]][17]),
          list(temp[[2]][1:2], temp[[2]][3], temp[[2]][4:9], temp[[2]][10:11], temp[[2]][12:14], temp[[2]][15], temp[[2]][16:17]),
          list(temp[[2]][1:2], temp[[2]][3], temp[[2]][4:11], temp[[2]][12:14], temp[[2]][15], temp[[2]][16:17]))
H2_s = lapply(H2, function(x) get_s_tissue(x, 2))

H3 = list(list(temp[[3]][1:6], temp[[3]][7:10], temp[[3]][11:13], temp[[3]][14:15], temp[[3]][16],temp[[3]][17]),
          list(temp[[3]][1:2], temp[[3]][3:4], temp[[3]][5:6], temp[[3]][7:8], temp[[3]][9], temp[[3]][10], temp[[3]][11:13], temp[[3]][14], temp[[3]][15:17]),
          list(temp[[3]][1:4], temp[[3]][5:6], temp[[3]][7:9], temp[[3]][10], temp[[3]][11:13], temp[[3]][14:15], temp[[3]][16], temp[[3]][17]),
          list(temp[[3]][1:6], temp[[3]][7:9], temp[[3]][10], temp[[3]][11:13], temp[[3]][14:15],  temp[[3]][16], temp[[3]][17]))
H3_s = lapply(H3, function(x) get_s_tissue(x, 3))

H4 = list(list(temp[[4]][1:6], temp[[4]][7:10], temp[[4]][11:13], temp[[4]][14:16]),
          list(temp[[4]][1], temp[[4]][2:5], temp[[4]][6], temp[[4]][7:9], temp[[4]][10], temp[[4]][11], temp[[4]][12:13], temp[[4]][14], temp[[4]][15], temp[[4]][16]),
          list(temp[[4]][1:6], temp[[4]][7:8], temp[[4]][9], temp[[4]][10], temp[[4]][11], temp[[4]][12], temp[[4]][13], temp[[4]][14:16]),
          list(temp[[4]][1:6], temp[[4]][7:8], temp[[4]][9], temp[[4]][10], temp[[4]][11:13], temp[[4]][14:16]))
H4_s = lapply(H4, function(x) get_s_tissue(x, 4))

H5 = list(list(temp[[5]][1], temp[[5]][2], temp[[5]][3:11], temp[[5]][12:17]),
          list(temp[[5]][1:3], temp[[5]][4:6], temp[[5]][7], temp[[5]][8:9], temp[[5]][10:11], temp[[5]][12:13], temp[[5]][14], temp[[5]][15], temp[[5]][16], temp[[5]][17]),
          list(temp[[5]][1], temp[[5]][2], temp[[5]][3:9], temp[[5]][10:11], temp[[5]][12], temp[[5]][13:17]),
          list(temp[[5]][1], temp[[5]][2], temp[[5]][3:11], temp[[5]][12], temp[[5]][13:17]))
H5_s = lapply(H5, function(x) get_s_tissue(x, 5))

H6 = list(list(temp[[6]][1:4], temp[[6]][5:8], temp[[6]][9:10], temp[[6]][11], temp[[6]][12:13],temp[[6]][14], temp[[6]][15:16],temp[[6]][17]),
          list(temp[[6]][1], temp[[6]][2:4], temp[[6]][5], temp[[6]][6], temp[[6]][7], temp[[6]][8], temp[[6]][9:11], temp[[6]][12:13],temp[[6]][14], temp[[6]][15:16],temp[[6]][17]),
          list(temp[[6]][1:4], temp[[6]][5:7], temp[[6]][8], temp[[6]][9:10], temp[[6]][11], temp[[6]][12:13],temp[[6]][14], temp[[6]][15:16],temp[[6]][17]),
          list(temp[[6]][1:4], temp[[6]][5:7], temp[[6]][8], temp[[6]][9:10], temp[[6]][11], temp[[6]][12:13],temp[[6]][14], temp[[6]][15:16],temp[[6]][17]))
H6_s = lapply(H6, function(x) get_s_tissue(x, 6))

Sil=as.data.frame(cbind(H1_s, H2_s, H3_s, H4_s, H5_s, H6_s))
Sil1=as.data.frame(cbind(H1_s, H2_s, H3_s, H4_s, H5_s, H6_s))

get_s_tissue <- function(clusters_heatmap, histone_type){
  dist_matrix = all_distances[[histone_type]]
  all_s = c()
  for (tissue1 in temp[[histone_type]]){
    print(tissue1)
    b_tissues = c()
    for (cluster in clusters_heatmap){
      if (tissue1 %in% cluster){
        if (length(cluster) == 1) a_tissue = 0 
        else a_tissue = mean(sapply(cluster, function(tissue2) return(dist_matrix[tissue1, tissue2])), na.rm=TRUE)
      } else {
        b_tissue = mean(sapply(cluster, function(tissue2) return(dist_matrix[tissue1, tissue2])))
        b_tissues = c(b_tissues, b_tissue)
      }
    }
    s = (min(b_tissues) - a_tissue)/ max(a_tissue, min(b_tissues))
    all_s = c(all_s, s)
  }
  names(all_s) = temp[[histone_type]]
  mean_s = mean(all_s[all_s != 1])
  return(mean_s)
}
clusters = list(epigenomes_annot$official_name[c(16,2,1,15,13,14,18)],
     epigenomes_annot$official_name[c(19)],
     epigenomes_annot$official_name[c(7,17,3,4)],
     epigenomes_annot$official_name[c(5,6,8,9,10,11,12)])  
temp = get_s_tissue(clusters, 1)
#======Rand index======
library("fossil")


all_distances = vector("list")
for (i in 1:length(all_sims)){
  all_distances[[i]] = 1-all_sims[[i]]
}

reorderfun = function(d, w) reorder(d, w)

d = dist(all_distances[[2]])
h = hclust(d, "complete")
d_de = as.dendrogram(h)
d_de = reorderfun(d_de, h$order)
d_de = as.hclust(d_de)

plot(d_de)
d_de_r = rect.hclust(d_de, k=2)

Hs = list( #H1
  list(c(rep(1, 7), rep(2, 1), rep(1, 4), rep(2, 7)),
       # c(1,2,2,1,1,2,1,3,2,1,2,2,3,3,1,2,4,2,3),
       c(2,1,1,2,2,1,2,3,1,2,1,1,3,3,2,1,4,1,3),
       c(rep(1, 7), 2, rep(1, 2), rep(3, 2), 2, 4, rep(2, 5)),
       c(rep(1, 7), 2, rep(1, 4), 2, 3,rep(2, 5))
  ), #H2
  list(c(rep(1, 2), 2, rep(1, 8), rep(2, 6)),
       c(1,2,3,2,1,2,1,2,2,1,1,3,4,1,3,2,3),
       # c(rep(1, 2), 2, rep(1,6), rep(3, 2), rep(2,3), 4, rep(2,2)),
       c(rep(2, 2), 3, rep(2,6), rep(4, 2), rep(3,3), 4, rep(3,2)),
       c(rep(1, 2), 2, rep(1,8), rep(2,3), 3, 2,2)
  ), #H3
  list(c(rep(1, 6), rep(2, 4), rep(1, 3), rep(2, 2), 1, 2),
       # c(rep(1, 2), rep(2, 2), rep(1, 2), rep(3,2), 2 , 3, rep(2,3), 3, rep(1,3)),       
       c(rep(3, 2), rep(1, 2), rep(3, 2), rep(2,2), 1 , 2, rep(1,3), 2, rep(3,3)),
       c(rep(1, 4), rep(2, 2), rep(3, 3), 4, rep(1, 3), rep(3, 2), 1, 3),
       c(rep(1, 6), rep(2, 3), 3, rep(1, 3), rep(2, 2), 1, 2)
  ), #H4
  list(c(rep(1, 6), rep(2, 4), rep(1, 3), rep(2, 3)),
       # c(1, rep(2, 4), 1, rep(3, 3), 1, 2, rep(1, 2), 3, 1, 2),
       c(2, rep(1, 4), 2, rep(3, 3), 2, 1, rep(2, 2), 3, 2, 1),
       c(rep(1, 6), rep(2, 2), 3, 2, 1, 4, 1, rep(2, 3)),
       c(rep(1, 6), rep(2, 2), 3, 2, rep(1, 3), rep(2, 3))
  ), #H5
  list(c(1, 2, rep(1, 9), rep(2, 6)),
       # c(rep(1, 3), rep(2, 3), 1, rep(2, 2), rep(1, 2), rep(3, 2), 1, 4, 2, 3),
       c(rep(2, 3), rep(1, 3), 2, rep(1, 2), rep(2, 2), rep(3, 2), 2, 4, 1, 3),
       # c(1, 2, rep(1, 7), rep(3, 2), 4, rep(2, 5)),
       c(1, 3, rep(1, 7), rep(2, 2), 4, rep(3, 5)),
       c(1, 2, rep(1,9), 3, rep(2, 5))
  ), #H6
  list(c(rep(1, 4), rep(2, 4), rep(1, 2), 2, rep(1, 2), 2,rep(1, 2), 2),
       # c(1, rep(2, 3), 3, 2, 1, 4, rep(1, 3), rep(2, 2), 4, rep(1, 2), 4),
       c(3, rep(1, 3), 2, 1, 3, 4, rep(3, 3), rep(1, 2), 4, rep(3, 2), 4),
       c(rep(1, 4), rep(2, 3), 3, rep(4, 2), 2, rep(1, 2), 2, rep(1, 2), 2),
       c(rep(1, 4), rep(2, 3), 3, rep(1, 2), 2, rep(1, 2), 2,rep(1, 2), 2)
  )
) 
# initial_clusters = clusters
clusters = list( #H1
  list(c(rep(1, 10), rep(2,9)), #Lifestage
       c(rep(1, 3), rep(2,7), rep(3,9)), #origin
       c(rep(1, 3), rep(2,7), rep(3,9)), #type
       c(rep(1, 10), rep(2,9)) #potency
  ), #H2
  list(c(rep(1, 2), rep(2,15)),
       c(rep(1, 2), rep(2,8), rep(3,7)),
       c(rep(1, 2), rep(2,8), rep(3,7)),
       c(rep(1, 2), rep(2,15))
  ), #H3
  list(c(rep(1, 4), rep(2,13)),
       c(rep(1, 4), rep(2,6), rep(3,7)), #Only 3 cats in Origin H3K36me3
       c(rep(1, 4), rep(2,6), rep(3,7)),
       c(rep(1, 4), rep(2,13))
  ), #H4
  list(c(rep(1, 6), rep(2,10)),
       c(rep(1, 6), rep(2,3), rep(3,7)), #Only 3 cats in Origin H3K4me1
       c(rep(1, 6), rep(2,3), rep(3,7)),
       c(rep(1, 6), rep(2,10))
  ), #H5
  list(c(rep(1, 9), rep(2,8)),
       c(rep(1, 9), rep(2,2), rep(3,6)),
       c(rep(1, 9), rep(2,2), rep(3,6)),
       c(rep(1, 9), rep(2,8))
  ), #H6
  list(c(rep(1, 4), rep(2,13)),
       c(rep(1, 4), rep(2,3), rep(3,10)),
       c(rep(1, 4), rep(2,3), rep(3,10)),
       c(rep(1, 4), rep(2,13))
  )
)
get_rand_adj <- function(Hs, clusters){
  rand.adj_res = list()
  for (i in 1:6){
    H_rand.adj = c()
    for (j in 1:4){
      H_rand.adj = c(H_rand.adj, round(adj.rand.index(Hs[[i]][[j]], clusters[[i]][[j]]), 3))
    }
    rand.adj_res[[i]] = H_rand.adj
  }
  names(rand.adj_res) = histone_type_list
  rand.adj_res = as.data.frame(rand.adj_res)
  rownames(rand.adj_res) = c("Life stage", "Origin", "Type", "Potency")
  rand.adj_res = rand.adj_res[c( "Potency", "Type", "Origin", "Life stage"),]
  return(as.data.frame(rand.adj_res))
} 

rand_adj = get_rand_adj(Hs, clusters)
fwrite(rand_adj, "rand_index.csv", sep='\t', quote=FALSE, dec = ',', row.names = TRUE, col.names = TRUE)
#------------New clustering--------------------
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

head(all_res_list.pearcor_padj_sig[[1]])
all_genes_H1 = Reduce(union, all_res_list.pearcor_padj_sig[[1]])
gene_cluster_H1 = lapply(all_genes_H1, function(y) lapply(all_res_list.pearcor_padj_sig[[1]], function(x)  y %in% x))
names(gene_cluster_H1) = all_genes_H1
gene_cluster_col = sapply(gene_cluster_H1, function(x) {
  all_tissues = names(x[x == TRUE])
  return(paste(all_tissues, collapse = ','))})
a = gene_cluster_col[1]
names(ZZ)

head(gene_cluster_col)
c1 = gene_cluster_col[1]

all_genes_clusters = lapply(c(gene_cluster_col[1:5], c1), function(gene_cluster) {
  gene_name = names(gene_cluster)
  print("======================================")
  all_tissues = Reduce(union, sapply(str_split(gene_cluster, ',')[[1]], function(x) str_split(x, '_')))
  print(all_tissues)
  adj_mat = get_adj_mat(gene_cluster)
  all_tissues_combi = unlist(Map(combn, list(all_tissues), seq(3, length(all_tissues)), simplify = FALSE), recursive=FALSE)
  all_tissues_combi = all_tissues_combi[order(sapply(all_tissues_combi, length), decreasing=T)]
  gene_cluster_list = check_cluster(all_tissues_combi)
  return(gene_cluster_list)
})
all_genes_clusters
temp = unlist(gene_cluster_col[1:5])
namestemp[2]
temp_list = c(gene_cluster_col[[1]], c1)
all_results = vector("list")
for (h in 1:length(temp_list)){
  gene_cluster = temp_list[h]
  all_tissues = Reduce(union, sapply(str_split(gene_cluster, ',')[[1]], function(x) str_split(x, '_')))
  print(all_tissues)
  adj_mat = get_adj_mat(gene_cluster)
  all_tissues_combi = unlist(Map(combn, list(all_tissues), seq(3, length(all_tissues)), simplify = FALSE), recursive=FALSE)
  all_tissues_combi = all_tissues_combi[order(sapply(all_tissues_combi, length), decreasing=T)]
  gene_cluster_list = check_cluster(all_tissues_combi)
  all_results[[h]] = gene_cluster_list
}
all_results

grep("AK9", names(all_genewise_cluster[[1]]))
all_genewise_cluster[[1]]["AOX1"]
names(all_genewise_cluster[[1]])[6652]

c1 = "a_d,a_c,c_d,b_d,b_c,d_e,a_e,c_e"
all_tissues = c("a","b","c","d","e")
adj_mat = get_adj_mat(c1, all_tissues)
Map(combn, list(c("a","b","c","d","e")), seq(3, length(c("a","b","c","d","e"))), simplify = FALSE)
temp = unlist(Map(combn, list(c("a","b","c","d","e")), seq(3, length(c("a","b","c","d","e"))), simplify = FALSE), recursive=FALSE)
all_tissues_combi = temp[order(sapply(temp, length), decreasing=T)]
temp = check_cluster(all_tissues_combi, adj_mat)
temp

