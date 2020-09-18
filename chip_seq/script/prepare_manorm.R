library(data.table)
library(dplyr)

setwd("/Users/dhthutrang/Documents/BIOINFO/Episplicing/ENCODE_episplicing/")
peak_metadata = fread("/Users/dhthutrang/Documents/BIOINFO/Episplicing/ENCODE_episplicing/chip_seq/peak_files/metadata.tsv")
# peak_metadata[peak_metadata$V10 == "spleen" & peak_metadata$V29 == "2016-11-03", "V10"] = "spleen2016"

# submeta = peak_metadata[, c("V1", "V7", "V10", "V22")]
# colnames(submeta) = c("peak", "read", "tissue", "histone")

submeta = peak_metadata %>%
  group_by(`Biosample term name`, `Experiment target`) %>%
  mutate(histone_type = strsplit(`Experiment target`, '-')[[1]][1],
         tissue_type = gsub('[[:punct:]]|\ ', '', `Biosample term name`),
         file_names = paste(histone_type, '_', tissue_type, '.bed', sep=''), 
         peak = paste('/home/dhthutrang/ENCODE/chip_seq/peak_files/merged/', file_names, sep=''),
         read = paste('/home/dhthutrang/ENCODE/chip_seq/alignment_files/bed/', file_names, sep=''),
         read = gsub('trophoblast', 'trophoblastcell', read),
         combined_files = paste(peak, read, sep="\t")) %>%
  ungroup() %>%
  dplyr::select(combined_files, histone_type, tissue_type) %>%
  unique()

for (histone_type in unique(submeta$histone_type)){
  temp1 = submeta[submeta$histone_type == histone_type,]
  temp2 = as.data.frame(t(combn(temp1$combined_files, 2, FUN=c)))
  temp3 = as.data.frame(t(combn(temp1$tissue_type, 2, FUN=c)))
  temp4 = cbind(temp2, temp3)
  colnames(temp4) = c("file_1", "file_2", "tis_1", "tis_2")
  temp = temp4 %>%
    mutate(tis1 = dplyr::if_else(tis_1 < tis_2, true = tis_1, false = tis_2),
           tis2 = dplyr::if_else(tis_1 < tis_2, true = tis_2, false = tis_1),
           file1 = dplyr::if_else(tis_1 < tis_2, true = file_1, false = file_2),
           file2 = dplyr::if_else(tis_1 < tis_2, true = file_2, false = file_1)
    ) %>%
    dplyr::select(file1, file2, tis1, tis2)
  fwrite(temp, paste("chip_seq", histone_type, "all_pairs.txt", sep='/'), sep="\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
}





