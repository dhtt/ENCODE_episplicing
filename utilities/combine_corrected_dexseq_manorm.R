## ---------------------------
##
## Script name: DEXSeq_analysis.R
##
## Purpose of script: Combine the multiple correction results
##
## Author: Trang Do
##
## Copyright (c) Trang Do, 2018
## Email: dhttrang@bioinformatik.uni-saarland.de
##
## ---------------------------
##
## Notes: 
## - Description: The script gathers the multiple correction results from global DEXSeq analysis and Manorm2 analysis
## and generate a dataframe containing the information about corrected exons
##
## ---------------------------

library("data.table")
library("dplyr")
library("optparse")


# ==== PREPARATION ====
# Parse arguments for input/output settings
option_list <- list(
  make_option("--global_DEXSeq_result",
    type = "character",
    help = "path to the result from global DEXSeq analysis in RDS format",
    metavar = "character",
    default = "mRNA_seq/dexseqcount/correction/count/res/dxd.res.RDS"
  ),
  make_option("--manorm2_result",
    type = "character",
    help = "path to all MAnorm2 results in xlsx format specific for histone marks",
    default = "chip_seq", 
    metavar = "character"
  ),
  make_option("--general_analysis_results",
    type = "character",
    help = "path to the folder where the general results from the analysis are stored and shared between processes",
    metavar = "character",
    default = "general_analysis_results"
  )
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
global_DEXSeq_result <- opt$global_DEXSeq_result
manorm2_res_path <- opt$manorm2_result
general_analysis_results <- opt$general_analysis_results

# Set working directory and input paths
working_dir <- getwd()
histone_types <- c("H3K27ac", "H3K27me3", "H3K36me3", "H3K4me3", "H3K9me3")


# ==== GATHER RESULTS FROM THE MULTIPLE CORRECTION WITH MANORM2 ====
get_manorm2_result <- function(manorm2_result_paths) {
  #' Get genes/exons that have DHMs in the MAnorm2 analysis from all histone types
  #' 
  #' The ID exons with DEU follows the format gene:exon_number
  #' @param manorm2_result_paths list of paths to MAnorm2 results in xlsx format
  #' @return a dataframe with information on which DHMs type the exons have and if the exons are first exons 
  #'
  
  # Get MAnorm2 results for all histone types and insert those in different columns in a dataframe DHM_df with IDs
  res_df_list <- lapply(manorm2_result_paths, fread)
  DHM_list <- lapply(res_df_list[2:length(manorm2_result_paths)], function(x) {return(x[["V10"]])})
  DHM_df <- cbind(res_df_list[[1]], do.call(cbind, DHM_list))
  colnames(DHM_df) <- c("chr", "source", "type", "start", "end", "score", "strand", "quality", "metadata", histone_types)
  DHM_df <- DHM_df[DHM_df$type != "aggregate_gene", ]

  # Extract gene/exon ID and mark first exons from DHM_df
  metadata_df <- lapply(DHM_df$metadata, function(x) unlist(strsplit(x, '"'))[c(2, 6, 8)])
  metadata_df <- as.data.frame(do.call("rbind", metadata_df))
  colnames(metadata_df) <- c("gene_id", "exon_id", "first_exon")
  metadata_df$exon_id[grep("E", metadata_df$exon_id, invert = T)] <- paste("E", metadata_df$exon_id[grep("E", metadata_df$exon_id, invert = T)], sep = "")
  metadata_df$id <- paste(metadata_df$gene_id, metadata_df$exon_id, sep = ":")
  DHM_df <- cbind(DHM_df, metadata_df)

  # Tidy up manorm2 result
  DHM_df <- DHM_df %>%
    dplyr::select(-metadata) %>%
    dplyr::group_by(id) %>%
    mutate(first_exon = first_exon[!is.na(first_exon)]) %>%
    dplyr::ungroup() %>%
    mutate(
      H3K27ac = dplyr::if_else(H3K27ac == ".", FALSE, TRUE),
      H3K27me3 = dplyr::if_else(H3K27me3 == ".", FALSE, TRUE),
      H3K36me3 = dplyr::if_else(H3K36me3 == ".", FALSE, TRUE),
      H3K4me3 = dplyr::if_else(H3K4me3 == ".", FALSE, TRUE),
      H3K9me3 = dplyr::if_else(H3K9me3 == ".", FALSE, TRUE)
    )
  DHM_df$dhm <- Reduce("|", DHM_df[histone_types])

  return(DHM_df)
}

# ==== COMBINE MANORM2 AND GLOBAL DEXSEQ RESULTS ====
combine_corrected_result <- function(manorm2_result, global_dexseq_result, corrected_exon_filename){
  #' Combine the information from MAnorm2 and global DEXSeq multiple correction
  #' 
  #' Results from the multiple correction using MAnorm2 and global DEXSeq are combined together. For each exon, the 
  #' significant DEU/DHM are marked in separate columns. This combined result does not contain the first exon in any
  #' gene. Only exons with DEU and any DHM are retained. 
  #' @param manorm2_result path to folder containing MAnorm2 results
  #' @param global_dexseq_result global DEXSeq's result as a dataframe 
  #' @param sig_exon_filename name given to the list of significant exons to store in RDS format
  #'
  
  # Get genes and exons having significant DHM event after correction with MAnorm2
  manorm2_result_files <- list.files(path = manorm2_result, full.names = TRUE)
  manorm2_result <- get_manorm2_result(manorm2_result_files)
  
  # Get genes and exons having significant DEU event after correction with DEXSeq
  global_DEXseq_result <- data.frame(readRDS(global_dexseq_result))
  exons_after_correction <- unique(global_DEXseq_result$groupID[global_DEXseq_result$padj < 0.05 & !is.na(global_DEXseq_result$padj)])

  # Combined info from MAnorm2 and global DEXSeq by filtering genes without global DEXSeq significant DEU 
  corrected_result_df <- manorm2_result
  corrected_result_df$deu <- sapply(unlist(corrected_result_df$id), function(x) x %in% exons_after_correction)

  # Filter first exons
  corrected_result_df <- corrected_result_df[corrected_result_df$first_exon == FALSE, ]
  
  # Mark DEUs or DHMs as significant
  corrected_result_df[c(histone_types, "dhm", "deu")] <- lapply(c(histone_types, "dhm", "deu"), 
    function(x) as.factor(corrected_result_df[x][[1]])
    )

  saveRDS(corrected_result_df, corrected_exon_filename) 
}

combine_corrected_result(manorm2_result, global_DEXSeq_result, 
                         paste(general_analysis_results, "combined_df_final.RDS", sep = "/"))
