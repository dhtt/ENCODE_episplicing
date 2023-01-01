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
## - Description: Script to compare DEUs identified from the pairwise and combined DEXSeq analysis
## - Preceeding script: DEXSeq_analysis_combined.R
## - Succeeding script: combine_corrected_dexseq_manorm.R TODO
##
## ---------------------------

library(optparse, quietly = TRUE)

# ==== PREPARATION ====
# Parse arguments for input/output settings
option_list <- list(
  make_option(c("-freq", "--DEUsFreqTable"),
    type = "character",
    help = "path to DEXSeq DEU frequency table in RDS format. $ENCODE_EXP/script/DEXSeq_DEU_freq_table.RDS was used",
    metavar = "character"
  ),
  make_option(c("-res", "--CombinedDEXSeqResult"),
    type = "character",
    help = "path to the result from combined DEXSeq analysis in RDS format. $ENCODE_EXP/res_90perc/dxd.res_90.RDS was used",
    metavar = "character"
  )
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Set working directory and input paths
working_dir <- getwd()

# Read in frequency table from multiple_correction.R for stored DEUs information
freq_table <- readRDS(opt$DEUsFreqTable)

# Read in result from the DEXSeq multiple correction done by DEXSeq_analysis_combined.R
combined_res_df <- data.frame(readRDS(opt$CombinedDEXSeqResult))


# ==== COMPARING OUTCOMES FROM THE PAIRWISE AND COMBINED DEXSEQ ANALYSIS ====
get_sig_exon = function(DEU_table, df, sig_exon_file_name){
  #' Get significant DEU from global correction
  #' 
  #' Each ID in gene_exon_id follows the format "gene;exon_number".
  #' This method gathers only the gene names or exon number.
  #' @param DEU_table pairwise DEXSeq's result  
  #' @param df combined DEXSeq's result as a dataframe 
  #' @param sig_exon_file_name name given to the list of significant exons to store in RDS format
  #'
  
  # Get genes and exons used for correction 
  genes_before_correction <- unique(DEU_table$gene)
  exons_before_correction <- DEU_table$exon

  no_genes_before_correction <- length(genes_before_correction)
  no_exons_before_correction <-  length(exons_before_correction)

  # Get genes and exons having significant DEU event after correction
  exons_after_correction <- unique(df[df$padj < 0.05 & !is.na(df$padj), c('groupID', 'featureID')])
  genes_after_correction <- unique(exons_after_correction$groupID)
  exons_after_correction <- paste(exons_after_correction$groupID, exons_after_correction$featureID, sep=':')

  no_exons_after_correction <- length(exons_after_correction)
  no_genes_after_correction <- length(genes_after_correction)
  
  # Check the set of exons that are overlap or union between PAIRWISE AND COMBINED DEXSEQ ANALYSIS
  overlap_combined_pairwise <- intersect(exons_after_correction, exons_before_correction)
  union_combined_pairwise <- union(exons_after_correction, exons_before_correction)
  combined_vs_pairwise <- setdiff(exons_after_correction, exons_before_correction)
  pairwise_vs_combined <- setdiff(exons_before_correction, exons_after_correction)

  # Check the set of genes that are overlap or union between PAIRWISE AND COMBINED DEXSEQ ANALYSIS
  overlap_combined_pairwise_gene <- intersect(genes_after_correction, genes_before_correction)
  combined_vs_pairwise_gene <- setdiff(genes_after_correction, genes_before_correction)
  pairwise_vs_combined_gene <- setdiff(genes_before_correction, genes_after_correction)
  
  # Report 
  print(paste("Number of corrected genes/exon:", 
              no_genes_before_correction, '/', no_exons_before_correction))
  print(paste("Number of significant genes/exons from combined DEXSeq:", 
              no_genes_after_correction, '/', no_exons_after_correction))
  
  print(paste("DEU rate exon (Combined/Pairwise):", no_exons_after_correction/no_exons_before_correction))
  print(paste("DEU rate gene (Combined/Pairwise):", no_genes_after_correction/no_genes_before_correction))
  print(paste("Overlap genes: ", length(overlap_combined_pairwise)))
  print(paste("Genes significant in combined but not in pairwise:", length(combined_vs_pairwise)))
  print(paste("Genes significant in pairwise but not in combined:", length(pairwise_vs_combined)))
  print(paste("Overlap/Union Ratio for significant genes:", length(overlap_combined_pairwise)/length(union_combined_pairwise)))
  saveRDS(exons_after_correction, sig_exon_file_name)
  return(overlap_combined_pairwise)
}
sig_exons = get_sig_exon(freq_table, combined_res_df, "multiple_correction_sig_exon_75low.RDS") #TODO
