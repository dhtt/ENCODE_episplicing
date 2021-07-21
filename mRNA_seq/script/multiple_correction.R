library(data.table, quietly=TRUE)
library(dplyr, quietly=TRUE)

setwd("/home/dhthutrang/ENCODE/mRNA_seq/dexseqcount/res")
all_files = list.files(pattern = '*.csv')


all_DEU_genes = vector("list", length(all_files))
for (i in 1:length(all_files)){
  file = all_files[i]
  print(file)
  data = data.table(read.csv(file, sep='\t'))
  all_DEU_genes[[i]] = unique(data[abs(data[[7]]) >= 2.0 & data[[6]] <= 0.05, 'groupID'][[1]])
}
  
all_DEU_genes = unique(Reduce(union, all_DEU_genes))
all_DEU_genes = all_DEU_genes[grep('+', all_DEU_genes, fixed = TRUE, invert = TRUE)]
all_DEU_genes = all_DEU_genes[grep('_', all_DEU_genes, fixed = TRUE, invert = TRUE)]
all_DEU_genes = paste(all_DEU_genes, collapse = "';")
print(all_DEU_genes)
