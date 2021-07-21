library(data.table, quietly=TRUE)
library(dplyr, quietly=TRUE)

setwd("/home/dhthutrang/ENCODE/mRNA_seq/dexseqcount/res")
all_files = list.files(pattern = '*.csv')

all_DEU_genes = vector("list", length(all_files))
for (i in 1:length(all_files)){
  file = all_files[i]
  print(file)
  data = data.table(read.csv(file, sep='\t'))
  all_DEU_genes[[i]] = unique(data[abs(data[[7]]) >= 2.0 & data[[6]] <= 0.0001, 'groupID'][[1]])
}
  
all_DEU_genes = unique(Reduce(intersect, all_DEU_genes))
all_DEU_genes = all_DEU_genes[grep('+', all_DEU_genes, fixed = TRUE, invert = TRUE)]
all_DEU_genes = all_DEU_genes[grep('_', all_DEU_genes, fixed = TRUE, invert = TRUE)]

print('Before intersect: ')
print(length(all_DEU_genes))
ref_genes = readRDS('/home/dhthutrang/ENCODE/mRNA_seq/script/all_genes.2021.RDS') 
all_DEU_genes = intersect(all_DEU_genes, ref_genes)
print('After intersect: ')
print(length(all_DEU_genes))

print(paste(all_DEU_genes, collapse = "';"))
