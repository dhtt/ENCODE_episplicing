## ---------------------------
##
## Script name: select_gene_for_correction.R
##
## Purpose of script: Filter the genes used for DEXSeq multiple correction
##
## Author: Trang Do
##
## Copyright (c) Trang Do, 2018
## Email: dhttrang@bioinformatik.uni-saarland.de
##
## ---------------------------
##
## Notes: 
## - Description: Given a folder containing all DEXSeq pariwise comparisons in csv format,
##   this script collect all exons that have differential usage in any sample pair, plot 
##   the cumulative distribution for the occurrence of DEUs and filter out genes that have
##   DEUs occurring in a high number of comparisons
## - Preceeding script: execute_dexseq.sh
## - Succeeding script: DEXSeq_analysis_combined.R 
##
## ---------------------------

library(data.table, quietly = TRUE)
library(dplyr, quietly = TRUE)
library(ggplot2, quietly = TRUE)
library(reshape2, quietly = TRUE)
library(optparse, quietly = TRUE)
library(doMC)
doMC::registerDoMC(cores = 8)

get_name <- function(ID_list, option = "gene") {
  #' Get only gene names/ exon numbers from the list of ID_list
  #'
  #' Each ID in gene_exon_id follows the format "gene;exon_number".
  #' This method gathers only the gene names or exon number.
  #' @param gene_exon_list list of gene_exon_id
  #'
  
  if (option == "gene") {
    return(unlist(lapply(ID_list, function(x) strsplit(x, split = ";")[[1]][1])))
  } else if (option == "exon") {
    return(unlist(lapply(ID_list, function(x) strsplit(x, split = ";")[[1]][2])))
  } else {
    warning("Check option for get_name()")
  }
}
general_analysis_results <- "/home/dhthutrang/ENCODE/general_analysis_results" # TODO


# ==== PREPARATION ====
# Parse arguments for input/output settings
option_list <- list(
  make_option(c("-p", "--DEXSeqCountPath"),
    type = "character",
    help = "path to DEXseq count in csv format. $ENCODE_EXP/dexseqcount was used",
    metavar = "character"
  )
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Set working directory and input paths
DEXSeqresultpath <- paste(opt$DEXSeqCountPath, 'res', sep = "/")
all_files <- list.files(DEXSeqresultpath, pattern = "*.csv")


# ==== INSPECT DEU GENES/EXONS ABUNDANCE ====
# Store genes with the exons having differential usage in all_DEU_genes.RDS
# These exons are defined by LFC >= 2 and FDR-adjusted p-values <= 0.0001
all_DEU_genes <- foreach(i = 1:length(all_files), .combine = "c", .packages = c("dplyr")) %dopar% {
  data <- data.table(read.csv(all_files[i], sep = "\t"))
  sig_exon <- data[abs(data[[7]]) >= 2.0 & data[[6]] <= 0.0001, ]
  sig_exon$id <- paste(sig_exon$groupID, sig_exon$featureID, sep = ";")
  return(sig_exon$id)
}
saveRDS(all_DEU_genes, paste(general_analysis_results, "all_DEU_genes.RDS", sep = "/"))


# Store the gene name, exon number and the number of tissue pairs in which the DEU
# for this exon occurred. The ID are sorted by decreased frequency. Overlapped genes
# which contain "+" or "-" in their names are excluded
freq <- table(all_DEU_genes)
freq <- freq[order(freq, decreasing = TRUE)]
freq <- freq[grep("+", names(freq), fixed = TRUE, invert = TRUE)]
freq <- freq[grep("-", names(freq), fixed = TRUE, invert = TRUE)]
freq_gene <- get_name(ID_list = names(freq), option = "gene")
freq_exon <- get_name(ID_list = names(freq), option = "exon")
freq_table <- data.table(freq_exon = as.numeric(freq), gene = freq_gene, exon = freq_exon)
saveRDS(freq_table, paste(general_analysis_results, "DEXSeq_DEU_freq_table.RDS", sep = "/"))


# ==== SUBFIGURE S1A: CUMULATIVE DISTRIBUTION FOR TISSUE PAIRS NUMBER IN WHICH AN DEU OCCURRED ====
tiff(paste(general_analysis_results, "ecdf_dexseq_correction.tiff", sep = "/"),
  res = 300, units = "in", width = 3, height = 3
)
ggplot(freq_table, aes(freq_exon, legend = F)) +
  stat_ecdf(pad = FALSE, na.rm = TRUE, size = 0.5) +
  scale_color_manual(values = "blue") +
  geom_vline(xintercept = 25, linetype = "dashed", color = "red", size = 0.5) +
  geom_hline(yintercept = 0.9, color = "red", size = 0.5) +
  theme_light() +
  scale_x_continuous(breaks = seq(0, 100, 10)) +
  scale_y_continuous(breaks = seq(0, 1, 0.1)) +
  xlab("Number of tissue pairs") +
  ylab("Probability") +
  theme(legend.position = "none", panel.grid.minor = element_blank(), aspect.ratio = 1)
dev.off()


# ==== FILTER NON-UBIQUITOUS DEUS FOR GLOBAL DEXSEQ ANALYSIS (MULTIPLE TESTING CORRECTION) ====
# Retain genes with non-ubiquitous DEUs, which are DEUs occurring in less than 26 
# pairwise comparisons and belonging to 90% of DEUs with the least occurence.
perc <- 0.1
filtered_exons <- freq[1:(round(length(freq) * perc))]
filtered_genes <- unique(get_name(names(filtered_exons)))
retained_genes <- setdiff(freq_gene, filtered_genes)
print("Number of genes with non-ubiquitous DEUs")


# Store retained genes with non-ubiquitous DEUs in multiple_correction_genes.txt
sink(paste(general_analysis_results, "multiple_correction_genes.txt", sep = "/"))
print(paste(all_DEU_genes_sig, collapse = ","))
sink()


# Filter genes in count files out if they are in filtered_genes list
count_folder <- paste(opt$DEXSeqCountPath, 'count', sep = "/")
count_files <- list.files(count_folder, full.names = TRUE)
filtered_count_folder <- paste(opt$DEXSeqCountPath, 'correction', 'count', sep = "/")

for (file in count_files) {
  file_name <- strsplit(file, "/")[[1]]
  output_path <- paste(filtered_count_folder, file_name[length(file_name)], sep="/")
  
  count_file <- fread(file, sep = "\t", col.names = c('ID', 'count'))
  filtered_ind <- sapply(count_file$ID, function(x) strsplit(x, ":")[[1]][1]) %in% filtered_genes
  filtered_ind[grep("+", count_file$ID, fixed = T )] = TRUE
  filtered_ind[grep("-", count_file$ID, fixed = T )] = TRUE
  count_file <- count_file[!filtered_ind, ]
  fwrite(count_file, output_path, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
}
