library(data.table, quietly=TRUE)
library(dplyr, quietly=TRUE)

setwd("~/ENCODE/mRNA_seq/dexseqcount/normedcount")
gene_list = readRDS('~/ENCODE/mRNA_seq/script/gene_list.2021.RDS')
all_files = list.files(full.names = F)
all_dfs = lapply(all_files, fread)
names(all_dfs) = all_files
print(head(all_files))

df = Reduce(function(x, y) merge.data.frame(x, y, by = c('groupID', 'featureID')), all_dfs) 
df = df[df$groupID %in% gene_list, ]
print(head(df))

tissue_names = setdiff(unique(sapply(colnames(df), function(x) strsplit(x, split = '_')[[1]][1])), c('groupID', 'featureID'))
normed_count = lapply(tissue_names, function(x) return(rowMeans(df[, grep(x, colnames(df))]) ))
names(normed_count) = tissue_names
normed_count = as.data.frame(normed_count)

normed_count_scaled = as.data.frame(t(apply(normed_count, 1, function(x) (x/max(x)))))
normed_count_scaled[is.na(normed_count_scaled)] = 0
saveRDS(normed_count_scaled, "~/ENCODE/mRNA_seq/script/normed_count_scaled.RDS")

#Gene wise
TSG = normed_count_scaled
TSG$gene_id = df$groupID
TSG = TSG %>% 
  dplyr::group_by(gene_id) %>%
  dplyr::summarise_all(mean(., na.rm=T))
TSG_mat = TSG %>% dplyr::select(-gene_id)
TSG_mat = as.data.frame(apply(TSG_mat, 1, function(x)(sum(1 - x/max(x))/(length(x) -1))))
TSG_mat$gene_id = TSG$gene_id
TSG_gene = TSG_mat$gene_id[TSG_mat >= 0.75]
saveRDS(TSG_gene, '~/ENCODE/mRNA_seq/script/TSG_gene.RDS')

# #Exon wise
# TSI = as.data.frame(apply(normed_count_scaled, 1, function(x)(sum(1 - x/max(x))/(length(x) -1))))
# colnames(TSI) = 'value'
# TSI$gene_id = df$groupID
# TSI$exon_id = df$featureID
# head(TSI, 10)
# saveRDS(TSI, '~/ENCODE/mRNA_seq/script/TSI.RDS')
# 
# TSI_res = TSI %>% filter(value > 0.95) %>% select(gene_id) %>% unique() %>% unlist()
# saveRDS(TSI_res, '~/ENCODE/mRNA_seq/script/TSG.RDS')
# print(paste('Number of TS gene', length(unique(TSI$gene_id)) - length(TSI_res)))
# 
# 
