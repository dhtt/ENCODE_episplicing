library(data.table, quietly=TRUE)
library(dplyr, quietly=TRUE)

# setwd("/Users/dhthutrang/Documents/BIOINFO/Episplicing/ENCODE_episplicing/mRNA_seq/script/temp")
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
# all_DEU_genes = readRDS("all_DEU_genes.RDS")

freq = table(all_DEU_genes)
freq = freq[order(freq, decreasing=TRUE)]
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

freq = freq[names(freq)[grep('+', names(freq), fixed = TRUE, invert = TRUE)]]
freq = freq[names(freq)[grep('_', names(freq), fixed = TRUE, invert = TRUE)]]
filtered_exon = names(freq[1:round(length(freq)*0.2)])

all_DEU_genes = intersect(unique(all_DEU_genes), filtered_exon)
all_DEU_genes = unique(unlist(lapply(all_DEU_genes, function(x) strsplit(x, split = ';')[[1]][1])))
# all_DEU_genes = all_DEU_genes[grep('+', all_DEU_genes, fixed = TRUE, invert = TRUE)]
# all_DEU_genes = all_DEU_genes[grep('_', all_DEU_genes, fixed = TRUE, invert = TRUE)]

print('Before intersect: ')
print(length(all_DEU_genes))
ref_genes = readRDS('/home/dhthutrang/ENCODE/mRNA_seq/script/all_genes.2021.RDS') 
all_DEU_genes = intersect(all_DEU_genes, ref_genes)
print('After intersect: ')
print(length(all_DEU_genes))

sink("multiple_correction_genes.txt")
print(paste(all_DEU_genes, collapse = "';"))
sink()
