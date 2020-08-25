library(data.table)
library(dplyr)

setwd("/Users/dhthutrang/Documents/BIOINFO/Episplicing/ENCODE_episplicing/")
mrna_metadata = fread("/Users/dhthutrang/Documents/BIOINFO/Episplicing/ENCODE_episplicing/file_transfer/submetadata_mrna.tsv")

submeta = mrna_metadata[, c("V1", "V10")]
colnames(submeta) = c("ID", "tissue")
submeta = submeta %>%
  mutate(tissue = gsub(' ', '', tissue),
         name = paste(paste(tissue, ID, sep='_'), '_count.txt', sep='')
         ) %>%
  dplyr::select(name)
fwrite(submeta, "/Users/dhthutrang/Documents/BIOINFO/Episplicing/ENCODE_episplicing/mRNA_seq/epi_ids.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
