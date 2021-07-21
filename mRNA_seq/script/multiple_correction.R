library(data.table, quietly=TRUE)
library(dplyr, quietly=TRUE)

setwd("/home/dhthutrang/ENCODE/mRNA_seq/dexseqcount/res")
all_files = list.files(pattern = '*.csv')


all_DEU_genes = vector("list", length(all_files))
for (i in length(all_files)){
  file = all_files[i]
  data = data.table(read.csv(file, sep='\t'))
  all_DEU_genes[[i]] = unique(data[abs(data[[7]]) >= 2.0 & data[[6]] <= 0.05, 'groupID'][[1]])
}
  
all_DEU_genes = Reduce(union, all_DEU_genes)
all_DEU_genes = paste(all_DEU_genes, collapse = "';")
print(all_DEU_genes)
