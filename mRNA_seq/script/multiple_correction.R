library(data.table, quietly=TRUE)
library(dplyr, quietly=TRUE)
library(ggplot2, quietly=TRUE)
library(reshape2, quietly=TRUE)

setwd("/Users/trangdo/Documents/Episplicing/ENCODE_episplicing/mRNA_seq/script/")
setwd("/home/dhthutrang/ENCODE/mRNA_seq/dexseqcount/res")
all_files = list.files(pattern = '*.csv')

#all_DEU_genes = vector("list", length(all_files))
#for (i in 1:length(all_files)){
#  file = all_files[i]
#  print(file)
#  data = data.table(read.csv(file, sep='\t'))
#  sig_exon = data[abs(data[[7]]) >= 2.0 & data[[6]] <= 0.0001, ]
#  sig_exon = sig_exon %>%
#  mutate(id = paste(groupID, featureID, sep=';'))
#  all_DEU_genes[[i]] = sig_exon$id
#}
#all_DEU_genes = Reduce(c, all_DEU_genes)
#saveRDS(all_DEU_genes, "all_DEU_genes.RDS")
all_DEU_genes = readRDS("all_DEU_genes.RDS")

length(all_DEU_genes)
freq = table(all_DEU_genes)
freq = freq[order(freq, decreasing=TRUE)]
freq = freq[grep('+', names(freq), fixed = T, invert = T)] 
freq = freq[grep('_', names(freq), fixed = T, invert = T)] 
print("SUMMARY")
print(summary(freq))
print("HEAD")
print(head(freq))
# tiff("ecdf.pdf") 
# plot(ecdf(freq))
# dev.off()
#tiff("ecdf.pdf") 
#plot(hist(freq))
#dev.off()

freq_gene = unlist(lapply(names(freq), function(x) strsplit(x, split = ';')[[1]][1]))
freq_exon = unlist(lapply(names(freq), function(x) strsplit(x, split = ';')[[1]][2]))
freq_table = data.table(freq_exon = as.numeric(freq), gene = freq_gene, exon = freq_exon) %>%
  dplyr::group_by(gene) %>%
  dplyr::summarize(rank_avg = mean(freq_exon),
                   rank_sum = sum(freq_exon))
freq_table_melt = reshape2::melt(freq_table)

all_DEU_genes = freq_table[order(freq_table$rank_avg, decreasing = T), "gene"][[1]]
paste(all_DEU_genes[1 : round(length(all_DEU_genes)*0.25)+1], collapse = ',')
paste(all_DEU_genes[(round(length(all_DEU_genes)*0.25)+1) : length(all_DEU_genes)], collapse = ',')
length(all_DEU_genes[(round(length(all_DEU_genes)*0.25)+1) : length(all_DEU_genes)] )

tiff("ecdf_dexseq_correction.tiff", res=300, units = "in", width=5, height = 5)
ggplot(freq_table_melt[freq_table_melt$variable == "rank_avg", ], aes(value, col=variable, legend=F)) +
  stat_ecdf(pad = F, na.rm = T) +
  scale_color_manual(values = "blue") +
  geom_vline(xintercept = 12.5, color="grey", linetype="dashed") + 
  geom_hline(yintercept = 0.75, color='red') +
  theme_light() +
  theme(legend.position = "none")
dev.off()

top_25_high = freq_table %>% dplyr::top_frac(0.25, rank_avg) %>% select(gene) %>% unlist()
top_75_low = freq_table %>% dplyr::top_frac(-0.75, rank_avg) %>% select(gene) %>% unlist()
paste(top_75_low, collapse = ',')

sink("multiple_correction_genes_75low.txt")
print(paste(top_75_low, collapse = ","))
sink()
sink("multiple_correction_genes_25high.txt")
print(paste(top_25_high, collapse = ","))
sink()
#==================

all_genes_sig = unlist(lapply(all_DEU_genes, function(x) strsplit(x, split = ';')[[1]][1]))
all_genes_sig = all_genes_sig[grep('+', all_genes_sig, fixed = T, invert = TRUE)]
all_genes_sig = all_genes_sig[grep('_', all_genes_sig, fixed = T, invert = TRUE)]
print('All genes with at least 1 sig DEU: ')
print(length(all_genes_sig))
all_genes_sig_freq = table(all_genes_sig)
all_genes_sig_freq = all_genes_sig_freq[order(all_genes_sig_freq, decreasing=TRUE)]
round(length(all_genes_sig_freq)*0.1)
tail(all_genes_sig_freq)
paste(names(all_genes_sig_freq[1:round(length(all_genes_sig_freq)*0.1)]), collapse=',')
paste(names(all_genes_sig_freq[round(length(all_genes_sig_freq)*0.9): length(all_genes_sig_freq)]), collapse=',')

perc=0.1
freq = freq[names(freq)[grep('+', names(freq), fixed = TRUE, invert = TRUE)]]
freq = freq[names(freq)[grep('_', names(freq), fixed = TRUE, invert = TRUE)]]
filtered_exon_1 = names(freq[1 : round(length(freq)*perc)])
filtered_exon_9 = names(freq[round(length(freq)*perc)+1 : length(freq)])
filtered_exon_1_low = names(freq[freq == 1])

get_name = function(gene_exon_id){
  return(unlist(lapply(gene_exon_id, function(x) strsplit(x, split = ';')[[1]][1])))
}
filtered_exon_1 = get_name(filtered_exon_1)
filtered_exon_9 = get_name(filtered_exon_9)
filtered_exon_1_low = get_name(filtered_exon_1_low)

length(setdiff(filtered_exon_1, filtered_exon_1_low))
length(intersect(all_genes_sig, filtered_exon_1_low))
"TSPAN3" %in% filtered_exon_1

print(paste(setdiff(filtered_exon_1, filtered_exon_1_low), collapse = ','))
# all_DEU_genes = all_DEU_genes[grep('+', all_DEU_genes, fixed = TRUE, invert = TRUE)]
# all_DEU_genes = all_DEU_genes[grep('_', all_DEU_genes, fixed = TRUE, invert = TRUE)]

print('Before intersect: ')
print(length(filtered_exon_1))
ref_genes = readRDS('all_genes.2021.RDS') 
all_DEU_genes_sig = intersect(filtered_exon_1, ref_genes)
print('After intersect: ')
print(length(all_DEU_genes_sig))

sink("/home/dhthutrang/ENCODE/mRNA_seq/script/multiple_correction_genes.txt")
print(paste(all_DEU_genes_sig, collapse = ","))
sink()
length(ref_genes)
