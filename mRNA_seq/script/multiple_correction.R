library(data.table, quietly=TRUE)
library(dplyr, quietly=TRUE)

# setwd("/Users/dhthutrang/Documents/BIOINFO/Episplicing/ENCODE_episplicing/mRNA_seq/script/temp")
setwd("/home/dhthutrang/ENCODE/mRNA_seq/dexseqcount/res")
all_files = list.files(pattern = '*.csv')

#all_DEU_genes = vector("list", length(all_files))
#for (i in 1:length(all_files)){
  #file = all_files[i]
  #print(file)
  #data = data.table(read.csv(file, sep='\t'))
  #print(data)
  #sig_exon = data[abs(data[[7]]) >= 2.0 & data[[6]] <= 0.0001, ]
  #sig_exon = sig_exon %>% 
   # mutate(id = paste(groupID, featureID, sep=';'))
  #all_DEU_genes[[i]] = sig_exon$id
  #print(all_DEU_genes[[i]][1:10])
#}
#all_DEU_genes = Reduce(union, all_DEU_genes)
#saveRDS(all_DEU_genes, "all_DEU_genes.RDS")
all_DEU_genes = readRDS("all_DEU_genes.RDS")
freq = table(all_DEU_genes)
freq = freq[order(freq, decreasing=FALSE)]
print(freq[0:10])
all_DEU_genes = unique(all_DEU_genes)
all_DEU_genes = unique(unlist(lapply(all_DEU_genes, function(x) strsplit(x, split = ';')[[1]][1])))
all_DEU_genes = all_DEU_genes[grep('+', all_DEU_genes, fixed = TRUE, invert = TRUE)]
all_DEU_genes = all_DEU_genes[grep('_', all_DEU_genes, fixed = TRUE, invert = TRUE)]

print('Before intersect: ')
print(length(all_DEU_genes))
ref_genes = readRDS('/home/dhthutrang/ENCODE/mRNA_seq/script/all_genes.2021.RDS') 
all_DEU_genes = intersect(all_DEU_genes, ref_genes)
print('After intersect: ')
print(length(all_DEU_genes))

#print(paste(all_DEU_genes, collapse = "';"))
