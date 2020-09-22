library(data.table, quietly=TRUE)
library(tidyverse)
library(dplyr, quietly=TRUE)
library(plyr)
library(rtracklayer)
library(ggplot2)
library(reshape2)
setwd("/Users/dhthutrang/Documents/BIOINFO/Episplicing/ENCODE_episplicing/flank")

#---------------------------------------------
all_res_list.pearcor_p = readRDS("all_res_list.pearcor_p.RDS")
lapply(all_res_list.pearcor_p, dim)
all_res_list.pearcor_p[[1]]

get_tissue_spec_ref <- function(){
  RPKM = read.csv("~/Documents/BIOINFO/Episplicing/files/flank/57epigenomes.RPKM.pc", row.names=1, sep="")
  head(RPKM)
  epigenomes = c("E003", "E004", "E005", "E006", "E007", "E011","E012","E013", "E016", "E024", "E053","E054", "E065","E066","E071","E079","E094","E095", "E096", "E098", "E100","E105","E106", "E109","E113") #E022-E027
  RPKM = RPKM[, epigenomes]
  norm_RPKM = transpose(as.data.frame(apply(RPKM, 1, function(x) (x/max(x)))))
  head(norm_RPKM)
  TSI_RPKM = as.data.frame(apply(norm_RPKM, 1, function(x)((length(norm_RPKM)-sum(as.numeric(x)))/length(norm_RPKM))))
  TSI_RPKM$gene = rownames(RPKM)
  head(TSI_RPKM)
  ENS_list = TSI_RPKM[TSI_RPKM[[1]] >= 0.75, 2]
  
  ENS_gene_list = read.delim("~/Documents/BIOINFO/Episplicing/files/flank/Ensembl_v65.Gencode_v10.ENSG.gene_info.txt", header=FALSE)
  head(ENS_gene_list)
  colnames(ENS_gene_list) = c("ens", "chr", "start", "end", "strand", "feature", "symbol", "name")
  TSI_symbols = ENS_gene_list[ENS_gene_list$ens %in% ENS_list, "symbol"]
  return(TSI_symbols)
}
TSI_symbols = get_tissue_spec_ref()
epigenomes_names = c("H1 Cells", "Mesendoderm", "Trophoblast",  "Mesenchyma", "Neuronal Progenitor Cells",  "Endoderm", "Ectoderm", "Mesoderm", "HUES64 Cells", "ES-UCSF4", "Cortex-derived Neurospheres", "Ganglion Eminence-derived Neurospheres",
                     "Aorta", "Liver", "Brain Hippocampus Middle", "Esophagus", "Gastric", "Left Ventricle", "Lung", "Pancreas", "Psoas Muscle",   "Right Ventricle", "Sigmoid Colon", "Small Intestine", "Spleen")

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
all_res_list.pearcor_p_sig = get_all_res_list_sig(all_res_list.pearcor_p, "pearcor_p", p_sig=0.1)
lapply(all_res_list.pearcor_p_sig, length)
lapply(all_res_list.pearcor_p_sig, function(x) lapply(x, function(y) length(y)))
lapply(all_res_list.pearcor_p_sig, function(x) lapply(x, function(y) paste(y, collapse = ', ')))

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
tissue_type_list = get_tissue_spec_array(all_res_list.pearcor_p_sig)
histone_type_list = list("H3K4me1", "H3K4me3", "H3K9me3", "H3K27me3", "H3K36me3", "H3K27ac")
tissue_type_list[[2]]$H1
#START here
get_genes_for_tissue <- function(res_list){
  all_genes_joined = vector("list") #List of 6 histone, each has sig genes for 25 tissues
  for (i in 1:length(res_list)){ #for each histone type
    res = res_list[[i]]
    tissue_list = tissue_type_list[[i]]
    all_tissues = lapply(tissue_list, function(x) return(Reduce(union, res[x])))
    names(all_tissues) = names(tissue_list)
    all_genes_joined[[i]] = all_tissues
  }
  return(all_genes_joined)
}
all_genes_joined = get_genes_for_tissue(all_res_list.pearcor_p_sig)
length(all_genes_joined)
lapply(all_genes_joined, length)
lapply(all_genes_joined[[5]], length)
lapply(all_genes_joined, function(x) lapply(x, function(y) length(y)))

merge_all <- function(x){ return(merge(x, by=0, all = T))}
get_all_len_before <- function(all_genes_joined){
  all_len_before = lapply(all_genes_joined, function(x) lapply(x, function(y) length(y)))
  all_len_before = lapply(all_len_before, function(x) as.data.frame(cbind(names(x), do.call(rbind, x))))
  all_len_before = all_len_before %>% reduce(full_join, by = "V1")
  colnames(all_len_before) = c('Tissue', histone_type_list)
  return(all_len_before)
}
all_len_before = get_all_len_before(all_genes_joined)

#============THIS PART HAVE NOT BEEN DONE============
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
  } #List of 6 histone, each has FILTERED sig genes for 25 tissues
  return(all_genes_joined)
}
all_genes_joined_filtered = get_genes_for_tissue_filtered(all_genes_joined)
length(all_genes_joined) #gene names
lapply(all_genes_joined, length)
lapply(all_genes_joined[[7]], length) 

get_all_len_after <- function(all_genes_joined){
  all_len_after = vector("list")
  for (i in 1:length(all_genes_joined)){
    len = lapply(all_genes_joined[[i]], length) 
    print(length(len))
    all_len_after[[i]] = len
  }
  all_len_after[[6]][c(10,11,12)] = NaN
  all_len_after = as.data.frame(do.call(cbind, all_len_after))
  rownames(all_len_after) = epigenomes_names
  colnames(all_len_after) = histone_type_list
  return(all_len_after)
}
all_len_after = get_all_len_after(all_genes_joined_filtered)
#============END here============

length(all_genes_joined_filtered[[1]])
paste(all_genes_joined_filtered[[4]][[12]], collapse = ', ')
temp = Reduce(intersect, all_genes_joined_filtered[[7]])
paste(temp, collapse = ', ')
temp = union(all_genes_joined_filtered[[7]][[12]], all_genes_joined_filtered[[2]][[12]])

# ----- Heatmap of tissue-spec genes -----
length(intersect(all_res_list.pearcor_p_sig_joinedtissue[[1]][[1]], all_res_list.pearcor_p_sig_joinedtissue[[1]][[2]]))
Reduce(intersect, all_res_list.pearcor_p_sig_joinedtissue[[5]])

all_distances = vector("list")
for (i in 1:length(all_genes_joined_filtered)){
  temp = all_genes_joined_filtered[[i]]
  distances = sapply(temp, function(x) sapply(temp, function(y) length(intersect(x, y))))
  diag(distances) = NA
  distances = distances/max(distances, na.rm = TRUE)
  rownames(distances) = epigenomes_names
  colnames(distances) = epigenomes_names
  if (i == 6){
    distances = distances[-c(10, 11, 12), -c(10, 11, 12)]
  }
  all_distances[[i]] = distances
}

library(NMF)
epigenomes_names
epigenomes_potency = c("Pluripotent", rep("Multipotent", 7), rep("Pluripotent", 2), 
                       rep("Multipotent", 2), rep("Differentiated", 13))
epigenomes_type = c("Primary culture", rep("ES cells-dev", 7), rep("Primary culture", 4), 
                    rep("Primary tissue", 13))
epigenomes_origin = c(rep("H1 cells", 5), rep("HUES64 cells", 4), 
                      "ES-UCSF4 cells", rep("Other", 15))
epigenomes_anatomy = c("ES cells", rep("ES cells-dev", 7), rep("ES cells", 2),
                       rep("Neurospheres", 2), "Heart", "Other", "Brain", 
                       "Digestive", "Digestive", "Heart", "Other", "Other",
                       "Muscle", "Heart", "Digestive", "Digestive", "Other")
epigenomes_annot = data.frame(epigenomes_potency, 
                              epigenomes_type,
                              epigenomes_origin, 
                              epigenomes_anatomy,
                              row.names = epigenomes_names)
head(epigenomes_annot)
lapply(epigenomes_annot, unique)
colnames(epigenomes_annot) = c("Potency", 
                               "Sample type",
                               "Origin",
                               "Anatomical location")
epigenomes_colors = list(c("lightgrey", "hotpink", "black"),
                         c("hotpink","black", "lightgrey"),
                         c("red", "lightgrey", "hotpink", "black"),
                         c("#FEF438", "#FC52B1","#FBC348", "#D23CC3", "#01075E", "#8A25EF", "#3B1FDB", "#0618B4")
)
# c("hotpink", "PaleTurquoise", "MediumPurple"), #old 
# c("pink", "plum","wheat", "salmon", "slateblue", "darkseagreen", "gold", "lightblue") #old
# )

tiff("pairwise_distance.tiff", width = 30, height = 40, units="in", res=100)
par(mfrow = c(4,2))
for (i in 1:7){
  if (i != 6){
    aheatmap(all_distances[[i]], Rowv = TRUE, Colv = TRUE, scale="none", 
             main=paste(histone_type_list[[i]]), annCol = epigenomes_annot, annRow = epigenomes_annot,
             annColors = epigenomes_colors, breaks = breaks, color = color )
  }
  else {
    aheatmap(all_distances[[i]], Rowv = TRUE, Colv = TRUE, scale="none", 
             main=paste(histone_type_list[[i]]), 
             annCol = epigenomes_annot[-c(10,11,12),], annRow = epigenomes_annot[-c(10,11,12),],
             annColors = list(epigenomes_colors[[1]], 
                              # epigenomes_colors[[2]][-c(1)]),
                              epigenomes_colors[[2]],#4 types
                              epigenomes_colors[[3]][-c(1)],#4 types
                              epigenomes_colors[[4]][-c(7)]), #4 types
             breaks = breaks, color = color)
  }
}
dev.off()

#separate images
folder = "/Users/dhthutrang/Documents/BIOINFO/Episplicing/episplicing/utilities/figures"
breaks = c(seq(0, 0.2, 0.2), seq(0.21, 1, 0.01))
color = '-RdYlBu2:82'
tiff(paste(paste(folder, histone_type_list[[1]], sep='/'), "tiff", sep='.'), width = 15, height = 12, units="in", res=100)
aheatmap(all_distances[[1]], Rowv = TRUE, Colv = TRUE, scale="none", 
         main = paste(histone_type_list[[1]]), annCol = epigenomes_annot, annRow = epigenomes_annot,
         annColors = epigenomes_colors, breaks = breaks, color = color, 
         legend = TRUE, annLegend = TRUE, labRow = NA, cellwidth = 25, cellheight = 15, fontsize = 20)
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
# ----- Annotation of common genes -----
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

cluster = list(list(1, c(2, 8, 7, 3), c(23, 22, 25, 17, 13, 15, 16, 24), c(6, 5, 4), c(21, 19, 20, 18)),
               list(2, c(18, 23,17,25), c(5, 2, 4), c(24, 20, 19, 16, 22, 21)),
               list(3, c(23, 13, 17, 15, 22, 21, 5, 24), c(16, 20, 25, 19), c(4, 2, 6), c(1, 12, 10, 11, 9)),
               list(4, c(1, 12, 10, 11), c(16, 19, 17, 23), c(3, 4, 6, 2), c(22, 21, 13, 25, 24, 20,14, 18, 15)),
               list(5, c(22, 21, 23), c(8,3,7,6), c(25, 15, 19, 17), c(24, 13, 18)),
               list(6, c(14, 18, 25), c(5, 8, 7), c(17, 13, 24, 23, 16, 19), c(20, 22, 21)), 
               list(7, c(20, 23, 21, 18), c(2, 3, 5, 4), c(22, 15, 17, 16), c(14, 24, 25, 19, 13)))
get_heatmap_subset <- function(reduce_method = "union"){
  genes_all_clusters = vector("list", length(cluster))
  for (j in 1:length(cluster)){
    histone_type = cluster[[j]][[1]]
    print(histone_type_list[histone_type])
    
    gene_each_cluster = vector("list", length(cluster[[j]]) - 1)
    for (i in 2:length(cluster[[j]])){
      group = cluster[[j]][[i]]
      group_name = paste(epigenomes_names[group], collapse = ', ')
      
      if (reduce_method == "union"){
        subset = Reduce(union, all_genes_joined_filtered[[histone_type]][group])
      }
      else if (reduce_method == "intersect"){
        subset = Reduce(intersect, all_genes_joined_filtered[[histone_type]][group])
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
genes_all_clusters = get_heatmap_subset("intersect")
head(genes_all_clusters[[1]])
lapply(genes_all_clusters, function(y) lapply(y, function(x) paste(x, collapse = ', ')))


lapply(cluster, function(x) paste(length(Reduce(union, all_genes_joined_filtered[[1]][x]))))
lapply(cluster, function(x) paste(Reduce(intersect, all_genes_joined[[1]][x]), collapse = ', '))

temp = Reduce(intersect, all_genes_joined_filtered[[5]][c(8,3,7,6)])
paste(temp, collapse = ", ")

Reduce(intersect, all_res_list.pearcor_p_sig_joinedtissue[[5]])

# -----3-----
all_len_after.his = all_len_after[,1:6]
all_len_after.his = as.data.frame(lapply(all_len_after.his, function(x) as.numeric(as.character(x))))
rownames(all_len_after.his) = epigenomes_names
colnames(all_len_after.his) = histone_type_list[1:6]
head(all_len_after.his)

n_tissues = 25
all_res_list_sig_joined = vector("list")
for (i in 1:n_tissues){
  all_res_sig = vector("list")
  for (j in 1:length(all_res_list.pearcor_p_sig_joinedtissue)){
    print( paste(j, i, sep = ", "))
    all_res_sig[[j]] = all_res_list.pearcor_p_sig_joinedtissue[[j]][[i]]
  }
  temp = Reduce(union, all_res_sig)
  all_res_list_sig_joined[[i]] = temp
}
length(all_res_list_sig_joined)
lapply(all_res_list_sig_joined, length)
lapply(all_res_list.pearcor_p_sig_joinedtissue[[6]][[10]], length)
all_res_list_sig_joined[[1]]

all_len_after.his$`Total genes with DHP (overlap)` = rowSums(all_len_after.his, na.rm=TRUE)
all_len_after.his$`Total genes with DHP (non-overlap)` = as.numeric(lapply(all_res_list_sig_joined, length))
sd_nan <- function(x){
  return(sd(x, na.rm = TRUE))
}
temp = all_len_after.his %>% 
  summarise_all(sd_nan) 

all_len_after.his$`Total genes with DMP (overlap)` = all_len_before$Methylation
all_len_after.his$`Total genes with DMP (non-overlap)` = all_len_after$Methylation

# -----FORMATTABLE TABLE-----
library(formattable)
head(all_len_before)
tissue_formatter <- formatter("span", 
                              style = x ~ style(
                                width = suffix(x, "px"),
                                font.weight = "bold", 
                                color = ifelse(x == "Total", "black", "gray")))
formattable(all_len_before,
            # align =c("l", "l", "c","c","c","c"), 
            list(
              `H3K4me1` = color_tile("white", "wheat"),
              `H3K4me3` = color_tile("white", "wheat"),
              `H3K9me3` = color_tile("white", "wheat"),
              `H3K27me3` = color_tile("white", "wheat"),
              `H3K36me3` = color_tile("white", "wheat"), 
              `H3K4me1` = color_tile("white", "wheat"),
              `H3K27ac` = color_tile("white", "wheat")
            )
)


