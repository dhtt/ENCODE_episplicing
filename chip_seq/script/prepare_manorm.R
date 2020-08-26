library(data.table)
library(dplyr)

setwd("/Users/dhthutrang/Documents/BIOINFO/Episplicing/ENCODE_episplicing/")
peak_metadata = fread("/Users/dhthutrang/Documents/BIOINFO/Episplicing/ENCODE_episplicing/file_transfer/submetadata_chipseq_peak.tsv")
peak_metadata[peak_metadata$V10 == "spleen" & peak_metadata$V29 == "2016-11-03", "V10"] = "spleen2016"

submeta = peak_metadata[, c("V1", "V7", "V10", "V22")]
colnames(submeta) = c("peak", "read", "tissue", "histone")
submeta = submeta %>%
  mutate(histone = sapply(histone, function(x) strsplit(x, split = "-")[[1]][1]),
         peak = paste('/home/dhthutrang/ENCODE/chip_seq/peak_files/', peak, '.bed', sep=''),
         read = paste('/home/dhthutrang/ENCODE/chip_seq/alignment_files/bed/', read, '.bed', sep=''),
         combined_files = paste(peak, read, sep="\t"),
         tissue = gsub(' ', '', tissue)) %>%
  dplyr::select(combined_files, tissue, histone)

for (histone_type in unique(submeta$histone)){
  temp1 = submeta[submeta$histone == histone_type,]
  temp2 = as.data.frame(t(combn(temp1$combined_files, 2, FUN=c)))
  temp3 = as.data.frame(t(combn(temp1$tissue, 2, FUN=c)))
  temp4 = cbind(temp2, temp3)
  fwrite(temp4, paste("chip_seq", histone_type, "all_pairs.txt", sep='/'), sep="\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
}