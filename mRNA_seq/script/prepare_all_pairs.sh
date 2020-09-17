library(data.table)
library(dplyr)
sam_folder = '/home/dhthutrang/ENCODE/mRNA_seq'
setwd("/home/dhthutrang/ENCODE/mRNA_seq/script")
mrna_metadata = fread("/home/dhthutrang/ENCODE/mRNA_seq/raw_data/metadata.tsv")

#submeta = mrna_metadata[, c("V1", "V10")]
#colnames(submeta) = c("ID", "tissue")

submeta = mrna_metadata %>%
  mutate(tissue = gsub('[[:punct:]]|\ ', '', `Biosample term name`),
	 filename = paste(sam_folder, '/', tissue, '_', `File accession`, '.sam', sep =''),
	 rename = paste('mv ', paste(sam_folder, '/', `File accession`, '.sam', sep=''), filename)
         ) %>%
  dplyr::select(rename)
fwrite(submeta, "/home/dhthutrang/ENCODE/mRNA_seq/script/renameSAMfile.sh", row.names = FALSE, col.names = FALSE, quote = FALSE)

tissue = mrna_metadata$`Biosample term name`
tissue = sapply(tissue, function(x) gsub('[[:punct:]]|\ ', '', x))
tissue = as.data.frame(unique(tissue))
fwrite(submeta, "/home/dhthutrang/ENCODE/mRNA_seq/script/epi_ids.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

