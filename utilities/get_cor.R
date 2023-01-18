## ---------------------------
##
## Script name: get_cor.R
##
## Purpose of script: -
##
## Author: Trang Do
##
## Copyright (c) Trang Do, 2018
## Email: dhttrang@bioinformatik.uni-saarland.de
##
## ---------------------------
##
## Notes:
## - Description: -
## - Preceeding script: -
## - Succeeding script: -
##
## ---------------------------

library("data.table", quietly = TRUE)
library("dplyr", quietly = TRUE)
library("boot", quietly = TRUE)
library("stats", quietly = TRUE)
library("parallel", quietly = TRUE)
library("tidyverse", quietly = TRUE)
library("optparse", quietly = TRUE)
library("doMC", quietly = TRUE)
doMC::registerDoMC(cores = 20)

# ==== PREPARATION ====
# Parse arguments for input/output settings
option_list <- list(
  make_option(c("-t", "--correlation_type"),
    type = "character",
    help = "use Pearson or Spearman correlation",
    default = "pearson",
    metavar = "character"
  ),
  make_option(c("-a", "--absolute_value"),
    type = "character",
    help = "convert M-values to absolute values",
    default = "FALSE"
  ),
  make_option(c("-d", "--dexseq_res_path"),
    type = "character",
    help = "path to DEXseq count results in csv format",
    default = "$ENCODE_EXP/dexseqcount/res",
    metavar = "character"
  ),
  make_option(c("-m", "--manorm_res_path"),
    type = "character",
    help = "path to MAnorm results in xlsx format",
    default = "$ENCODE_HIS", # ENCODE/chip_seq
    metavar = "character"
  ),
  make_option(c("-c", "--corrected_res_path"),
    type = "character",
    help = "path to multiple correction results in RDS format",
    default = "$ENCODE/utilities/combined_df_exon_90_final.RDS",
    metavar = "character"
  ),
  make_option(c("-i", "--all_pairs_DEU_path"),
    type = "character",
    help = "path to save the final dataframe containing DEU statistics",
    default = "/home/dhthutrang/ENCODE/flank/150123/all_pairs_DEU.RDS", # TODO
    metavar = "character"
  ),
  make_option(c("-j", "--all_pairs_DEU_flt_path"),
    type = "character",
    help = "path to save the filtered dataframe containing DEU statistics",
    default = "/home/dhthutrang/ENCODE/flank/150123/all_DEUs_corrected.RDS", # TODO
    metavar = "character"
  ),
  make_option(c("-k", "--all_pairs_DHM_path"),
    type = "character",
    help = "path to save the final dataframe containing DHM M-values",
    default = "/home/dhthutrang/ENCODE/flank/150123/all_pairs_DHM_list.RDS", # TODO
    metavar = "character"
  ),
  make_option(c("-l", "--all_pairs_DHM_flt_path"),
    type = "character",
    help = "path to save the filtered dataframe containing DHM M-values",
    default = "/home/dhthutrang/ENCODE/flank/150123/all_DHMs_corrected_list_manorm.RDS", # TODO
    metavar = "character"
  ),
  make_option(c("-p", "--p_values_path"),
    type = "character",
    help = "path to save p-values from correlation",
    default = "/home/dhthutrang/ENCODE/flank/150123/all_res_list.pearcor_p_90_manorm.RDS", # TODO
    metavar = "character"
  ),
  make_option(c("-r", "--r_values_path"),
    type = "character",
    help = "path to save r-values from correlation",
    default = "/home/dhthutrang/ENCODE/flank/150123/all_res_list.pearcor_r_90_manorm.RDS", # TODO
    metavar = "character"
  )
)
opt <- parse_args(OptionParser(option_list = option_list))

dexseq_res_path <- opt$dexseq_res_path
manorm_res_path <- opt$manorm_res_path
corrected_res_path <- opt$corrected_res_path
all_pairs_DEU_path <- opt$all_pairs_DEU_path
all_pairs_DEU_flt_path <- opt$all_pairs_DEU_flt_path
all_pairs_DHM_path <- opt$all_pairs_DHM_path
all_pairs_DHM_flt_path <- opt$all_pairs_DHM_flt_path
p_values_path <- opt$p_values_path
r_values_path <- opt$r_values_path

correlation_option_list <- list(opt$correlation_type, opt$absolute_value)
names(correlation_option_list) <- c("correlation_type", "absolute_value")
histone_type_list <- c("H3K27ac", "H3K27me3", "H3K36me3", "H3K4me3", "H3K9me3")


# ==== GENERIC FUNCTIONS ====
absmax <- function(x) {
  #' Get the value with maximum absolute value in an array/list x
  #'
  #' @param x a list of numbers
  #' @return value which has the maximum absolute value in x

  x[which.max(abs(x))]
}

get_colname <- function(filename_list, option = "DHM") {
  #' Get the list of sample name from file paths to use as column name for the combined result dataframe
  #'
  #' @param filename_list file paths that contain info for which 2 samples are being compared
  #' @param option either DEU or DHM, as the result dataframes of these have different shapes
  #' @return a list of sample pairs that can be used for aggregated result labeling
  #'

  name <- sapply(filename_list, function(x) strsplit(x, split = "/"))
  name <- sapply(name, function(x) x[length(x)][[1]])
  if (option == "DHM") {
    name <- sapply(name, function(x) strsplit(x, split = "\\.")[[1]][1])
  } else if (option == "DEU") {
    name <- sapply(name, function(x) strsplit(x, split = ".res.csv")[[1]][1])
  }
  column_name <- sapply(name, function(x) paste(sort(strsplit(x, split = "_")[[1]]), collapse = "_"))
  names(column_name) <- NULL
  return(column_name)
}

filter_genes <- function(result_df, filter_genes_path = "combined_df_exon.RDS", filter = "deu") {
  #' Filter out unsignificant DEU events from the dataframe containing all pairwise results using multiple correction
  #' outcome
  #'
  #' This function filters out first exon of any gene and genes that were excluded in the multiple comparisons
  #' correction step
  #' @param result_df a dataframe containing the results of pairwise comparisons for DEU/DHM. The first 2 columns should
  #' store gene ID and exon ID. From the 3rd column, the pairwise results
  #' are stored.
  #' @param filter_genes_path path to the dataframe that store the results from the global analysis for multiple testing
  #' correction. Columns are named after the features, according to which result_df would be filtered by. 
  #' filter_genes_path$type should contain only exonic_part, while first_exon contains both TRUE and FALSE values
  #' @param filter name of the feature the result_df should be filtered by, includes "deu" or any other histone type
  #' @return a data table with insigficant DEU/DHM set to 0 and without first exons
  #'

  # Get corrected genes for filtering from filter_genes_path file
  corrected_exons <- readRDS(filter_genes_path)
  corrected_genes <- unique(unlist(corrected_exons[corrected_exons["deu"] == TRUE, "gene_id"]))

  # Filter by corrected genes (10% of initial genes that have non-ubiquitous DEUs)
  result_df <- as.data.frame(result_df[result_df$gene_id %in% corrected_genes, ])
  corrected_df <- as.data.frame(corrected_exons[corrected_exons$gene_id %in% corrected_genes, ])

  # Order result_df and corrected_df so that their rownames overlap
  result_df <- result_df[order(result_df$gene_id), ]
  corrected_df <- corrected_df[order(corrected_df$gene_id), ]

  # Filter out first exons
  result_df <- result_df[as.logical(corrected_df$first_exon) == FALSE, ]
  corrected_df <- corrected_df[as.logical(corrected_df$first_exon) == FALSE, ]

  # Filter by DEU/DHM: If the DEU/DHM event was not significant as found in corrected_df,
  # the value in result_df is reset to 0
  final_filtered_df <- result_df
  if (filter == "deu") final_filtered_df[unlist(corrected_df[[filter]] == FALSE), 3:length(final_filtered_df)] <- 0

  return(as.data.table(final_filtered_df))
}

# ==== GET DEU FROM PAIRWISE DEXSEQ RESULTS ====
get_all_pairs_DEU <- function(DEU_list) {
  #' Parse results from all DEXSeq pairwise comparisons and store them into dataframe
  #'
  #' This function extracts the DEU statistics and FDR adjusted p-values from the DEXSeq pairwise comparison results,
  #' sets unsignificant DEUs' stats to 0, then saves the stats for all genes in each comparison into a column in the
  #' final dataframe.
  #' @param DEU_list list of paths to the results of pairwise comparisons
  #' @return a dataframe of DEU statistics with sample pairs as columns and gene/exons as rows
  #'

  colname_exp <- c("gene_id", "exon_id", get_colname(DEU_list, "DEU"))
  no_comparisons <- length(DEU_list)
  pair_DEU_list <- vector("list", no_comparisons)
  exp_id <- NA

  # Read DEXSeq pairwise comparison DEU statistics and FDR into a list pair_DEU_list
  for (i in seq(no_comparisons)) {
    print(paste("Pair: ", i, ": ", DEU_list[[i]], sep = ""))
    pair_DEU <- fread(DEU_list[[i]])
    exp_id <- pair_DEU[, 1:2]
    colnames(pair_DEU)[ncol(pair_DEU)] <- "LFC"
    pair_DEU_list[[i]] <- pair_DEU[, c("stat", "padj", "LFC")]
  }

  # For each result in pair_DEU_list, set stat=0 if DEU is not significant (FDR<0.05)
  pair_DEU_list <- lapply(
    pair_DEU_list,
    function(x) {
      x <- x %>%
        mutate(
          direction = LFC/abs(LFC),
          exp = dplyr::if_else(padj <= 0.05 & !is.na(padj), true = stat*direction, false = 0.0)
        ) %>%
        dplyr::select(exp)
    }
  )

  # Combine the DEU stats for the pairwise comparisons into 1 dataframe
  pair_DEU_list <- as.data.frame(pair_DEU_list)
  pair_DEU_list <- as.data.frame(cbind(exp_id, pair_DEU_list))
  pair_DEU_list <- pair_DEU_list[order(pair_DEU_list[[1]]), ]
  colnames(pair_DEU_list) <- colname_exp

  return(as.data.table(pair_DEU_list))
}

# Execute get_all_pairs_DEU()
print("AGGREGATING DEU RESULTS .....")
all_pairs_DEU <- list.files(dexseq_res_path, full.names = TRUE, pattern = "_res.csv")
all_pairs_DEU <- get_all_pairs_DEU(all_pairs_DEU)
saveRDS(all_pairs_DEU, all_pairs_DEU_path)

# Filter out unsignificant DEU events according to multiple correction results
print("FILTERING DEU RESULTS .....")
all_DEUs_corrected <- filter_genes(all_pairs_DEU, filter_genes_path = corrected_res_path, filter = "deu")
saveRDS(all_DEUs_corrected, all_pairs_DEU_flt_path)


# ==== GET DHM FROM PAIRWISE MANORM RESULTS ====
get_all_pairs_DHM <- function(DHM_list, his_type) {
  #' Parse results from all MAnorm pairwise comparisons for an histone type and store them into dataframe
  #'
  #' This function extracts the M-values and FDR adjusted p-values from the MAnorm pairwise comparison results, sets the
  #' M-value of unsignificant DHM event to 0, then saves M-values for all genes in each comparison into a column in the
  #' final dataframe. This function is called from get_all_pairs_DHM_list() for all investigated histone types.
  #' @param DHM_list list of paths to the results of pairwise comparisons
  #' @param his_type type of histone modification being processed
  #' @return a dataframe of DHM M-values with sample pairs as columns and gene/exons as rows
  #'

  pair_DHM_list <- vector("list", length(DHM_list))
  pair_DHM_list <- foreach(i = seq(length(DHM_list)), .combine = "cbind", .packages = c("dplyr")) %dopar% {
    # Read MAnorm pairwise comparison M-values and FDR of a comparison into pair_DHM
    print(paste("Pair: ", i, sep = ""))
    pair_DHM <- fread(DHM_list[[i]])
    id <- as.data.frame(do.call(rbind, lapply(pair_DHM$V9, function(x) strsplit(x, split = '"', fixed = T)[[1]][c(2, 6)])))
    colnames(id) <- c("gene", "exon")

    # Set M-values=0 if DHM is not significant (FDR<0.05)
    pair_DHM <- pair_DHM %>%
      dplyr::mutate(
        gene = id$gene, exon = id$exon, type = V3,
        p_val = as.numeric(as.character(V11)),
        m_val = dplyr::if_else(
          p_val <= 0.05,
          true = as.numeric(as.character(V10)),
          false = 0
        )
      ) %>%
      dplyr::select(gene, exon, m_val) %>%
      dplyr::group_by(gene, exon) %>%
      dplyr::summarise_all(absmax) %>%
      dplyr::na_if(., -Inf) %>%
      dplyr::na_if(., Inf) %>%
      dplyr::ungroup()

    # Save the M-values without gene/exon id
    if (i != 1){
      pair_DHM <- pair_DHM %>% dplyr::select(-gene, -exon)
    } 

    # Convert M-values if their absolute values should be used
    if (correlation_option_list$absolute_value == TRUE) {
      pair_DHM <- abs(pair_DHM)
    }

    # Store M-values of the current pairwise comparison into pair_DHM_list
    return(as.data.frame(pair_DHM))
  }
  print(dim(pair_DHM_list))
  return(pair_DHM_list)
}

get_all_pairs_DHM_list <- function(histone_type_list) {
  #' Perform results aggregation by calling get_all_pairs_DHM() for all histone types in the study
  #'
  #' @param histone_type_list list of histone types for which MAnorm results are available
  #' @return a list of dataframe of DHM M-values
  #'

  all_pairs_DHM_list <- vector("list", length(histone_type_list))

  # For each histone mark, parse, process and store the M-values from all pairwise comparisons.
  for (j in 1:length(histone_type_list)) {
    his <- histone_type_list[[j]]
    all_pairs_DHM <- list.files(paste(manorm_res_path, his, "flank/fl", sep = "/"), pattern = ".txt", full.names = TRUE)
    colname_his <- c("gene_id", "exon_id", get_colname(all_pairs_DHM, "DHM"))
    all_pairs_DHM.sig <- get_all_pairs_DHM(all_pairs_DHM, his)
    colnames(all_pairs_DHM.sig) <- colname_his

    all_pairs_DHM_list[[j]] <- as.data.table(all_pairs_DHM.sig)
  }
  return(all_pairs_DHM_list)
}

filter_all_his_list <- function(histone_df_list, histone_type_list, filter_genes_path) {
  #' Filter out unsignificant DHM events according to multiple correction results
  #'
  #' @param histone_df_list list of the aggregated results for all histone types from get_all_pairs_DHM_list()
  #' @param histone_type_list list of histone types for which MAnorm results are available
  #' @param filter_genes_path see filter_genes_path param from filter_genes()
  #' @return a list of dataframe with insigficant DEU/DHM set to 0 and without first exons
  #'

  all_filtered_df <- vector("list")

  # For each histone mark, filter out unsignificant DHM events using multiple correction results
  for (i in 1:length(histone_type_list)) {
    histone <- histone_type_list[i]
    his_df <- histone_df_list[[i]]
    all_filtered_df[[i]] <- filter_genes(result_df = his_df, filter_genes_path = filter_genes_path, filter = histone)
  }
  names(all_filtered_df) <- histone_type_list
  return(all_filtered_df)
}


# Execute get_all_pairs_DHM_list()
print("AGGREGATING DHM RESULTS .....")
all_pairs_DHM_list <- get_all_pairs_DHM_list(histone_type_list)
saveRDS(all_pairs_DHM_list, all_pairs_DHM_path)
all_pairs_DHM_list_ <- readRDS(all_pairs_DHM_path)
names(all_pairs_DHM_list_) <- histone_type_list

# Filter out unsignificant DEU events according to multiple correction results
print("FILTERING DHM RESULTS .....")
all_DHMs_corrected_list <- filter_all_his_list(histone_df_list =  all_pairs_DHM_list_, 
                                                 histone_type_list = histone_type_list,
                                                 filter_genes_path = corrected_res_path)
saveRDS(all_DHMs_corrected_list, all_pairs_DHM_flt_path)

print("===== FINISH AGGREGATING DATA =====")
# ===== PERFORM CORRELATION TEST =====
p_value_calculator <- function(r, nrow) {
  #' Calculate the p-values given a correlation coefficients. As reference, rcorr() from package Hmisc was used.
  #'
  #' @param r correlation coefficients
  #' @param nrow number of comparisons
  #' @return a p-value
  #'

  P <- r * sqrt(nrow - 2) / sqrt(1 - r * r)
  P <- 2 * pt(-abs(P), nrow - 2)
  return(P)
}

compute_cor_p <- function(DEU, DHM, cor_type = "pearson") {
  #' Compute Pearson's correlation p-value for a gene using its DEU and DHM values
  #'
  #' @param DEU DEXseq DEU statistics array
  #' @param DHM MAnorm M-values array
  #' @param cor_type Pearson or Spearman correlation
  #' @return a p-value
  #'

  # Only compute for gene having more than 1 unique DEU/DHM values
  if (length(unique(DEU)) > 1 & length(unique(DHM)) > 1) {
    p_val <- p_value_calculator(cor(DEU, DHM, method = cor_type), nrow = length(DEU))
    return(p_val)
  } else {
    return(NA)
  }
}

compute_cor_r <- function(DEU, DHM, cor_type = "pearson") {
  #' Compute Pearson's correlation coefficients for a gene using its DEU and DHM values
  #'
  #' @param DEU DEXseq DEU statistics array
  #' @param DHM MAnorm M-values array
  #' @param cor_type Pearson or Spearman correlation
  #' @return an R-value
  #'

  df <- as.data.frame(cbind(DEU, DHM))
  n_sep_point <- nrow(unique(df))

  # Only compute for gene having more than 1 unique DEU/DHM values and more than 2 unique data points (each data point
  # is defined by the pair of DEU DHM values)
  if (0 %in% apply(df, 1, unique)) n_sep_point <- n_sep_point - 1
  if (n_sep_point > 2 & length(unique(DEU)) > 1 & length(unique(DHM)) > 1) {
    r_val <- cor(DEU, DHM, method = cor_type)
    return(r_val)
  } else {
    return(NA)
  }
}

compute_cor <- function(all_pairs_DEU, all_pairs_DHM, compute_option = "p", all_genes, cor_type = "pearson") {
  #' Calculate the correlation coefficients and p-values for a gene using its DEU and one of the DHM values.
  #'
  #' @param all_pairs_DEU DEXseq DEU statistics array
  #' @param all_pairs_DHM MAnorm M-values array
  #' @param compute_option compute correlation coefficients "r" or p-values "p"
  #' @param all_genes list of genes for which correlation is computed. The order of genes in all_genes corresponds to
  #' the order of correlation results
  #' @param cor_type see cor_type param in compute_cor_p() and compute_cor_r()
  #' @return a dataframe containing p-values or r-values with sample pairs as columns and genes as rows
  #'

  # Since different DHM types have different availability of samples, the columns in DEU dataframe (all_pairs_DEU) is
  # filtered to contain only comparisons available for the currently processed DHM
  all_res_pair <- vector("list", ncol(all_pairs_DEU) - 2)
  subset_name <- colnames(all_pairs_DHM)
  colnames(all_pairs_DEU) <- gsub("trophoblastcell", "trophoblast", colnames(all_pairs_DEU))
  all_pairs_DEU_subset <- all_pairs_DEU[, ..subset_name]

  # Compute correlation coefficients or p-values for all sample pairs in parallel
  all_res_pair <- foreach(i = 1:(length(subset_name) - 2), .combine = "c", .packages = c("dplyr")) %dopar% {
    print(paste("Pair: ", i, sep = ""))
    exp <- all_pairs_DEU_subset[[i + 2]]
    his <- all_pairs_DHM[[i + 2]]
    data_table <- as.data.table(cbind(exp, his))

    if (compute_option == "p") {
      res_table <- data_table %>%
        dplyr::group_by(all_pairs_DEU$gene_id) %>%
        dplyr::summarise(res = compute_cor_p(exp, his, cor_type = cor_type)) %>%
        dplyr::select(res)
    } else if (compute_option == "r") {
      res_table <- data_table %>%
        dplyr::group_by(all_pairs_DEU$gene_id) %>%
        dplyr::summarise(res = compute_cor_r(exp, his, cor_type = cor_type)) %>%
        dplyr::select(res)
    }
  }
  all_res_pair <- as.data.table(all_res_pair)
  all_res_pair <- cbind(all_genes, all_res_pair)
  colnames(all_res_pair) <- c("gene_id", subset_name[3:length(subset_name)])

  return(as.data.frame(all_res_pair))
}

compute_cor_list <- function(all_pairs_DEU, all_pairs_DHM_list, compute_option = "p", cor_type = "pearson") {
  #' Perform compute_cor() for all types of DHM
  #'
  #' @param all_pairs_DEU DEXseq DEU statistics array
  #' @param all_pairs_DHM MAnorm M-values array
  #' @param compute_option compute correlation coefficients "r" or p-values "p"
  #' @param cor_type see cor_type param in compute_cor_p() and compute_cor_r()
  #' @return a list of dataframe containing p-values or r-valuesgenerated from compute_cor()
  #'

  all_res_list <- vector("list", length(histone_type_list) - 1)
  all_genes <- unique(all_pairs_DEU$gene_id)

  # Call compute_cor() for each histone type
  for (j in 1:length(histone_type_list)) {
    print(paste("Histone:", histone_type_list[[j]]))
    all_pairs_DHM <- all_pairs_DHM_list[[j]]
    if (compute_option == "p") {
      all_res_pair <- compute_cor(
        all_pairs_DEU = all_pairs_DEU,
        all_pairs_DHM = all_pairs_DHM,
        compute_option = "p",
        all_genes = all_genes,
        cor_type = cor_type
      )
    } else if (compute_option == "r") {
      all_res_pair <- compute_cor(
        all_pairs_DEU = all_pairs_DEU,
        all_pairs_DHM = all_pairs_DHM,
        compute_option = "r",
        all_genes = all_genes,
        cor_type = cor_type
      )
    }
    all_res_list[[j]] <- all_res_pair
  }
  return(all_res_list)
}

print("COMPUTING CORRELATION P-VALUES .....")
p_values_list <- compute_cor_list(
  all_pairs_DEU = all_DEUs_corrected,
  all_pairs_DHM_list = all_DHMs_corrected_list,
  compute_option = "p",
  cor_type = correlation_option_list$correlation_type
)
saveRDS(p_values_list, p_values_path)

print("COMPUTING CORRELATION R-VALUES .....")
r_values_list <- compute_cor_list(
  all_pairs_DEU = all_DEUs_corrected,
  all_pairs_DHM_list = all_DHMs_corrected_list,
  compute_option = "r",
  cor_type = correlation_option_list$correlation_type
)
saveRDS(r_values_list, r_values_path)
print("===== FINISH CORRELATION COMPUTATION =====")
