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

library(data.table, quietly = TRUE)
library(tidyverse, quietly = TRUE)
library(dplyr, quietly = TRUE)
library(plyr, quietly = TRUE)
library(rtracklayer, quietly = TRUE)
library(ggplot2, quietly = TRUE)
library(reshape2, quietly = TRUE)
library(NMF, quietly = TRUE)
library(formattable, quietly = TRUE)
library(viridis, quietly = TRUE)
library(ggplot2, quietly = TRUE)
library(ggraph, quietly = TRUE)
library(ggpubr, quietly = TRUE)
library(RColorBrewer, quietly = TRUE)
library(doMC, quietly = TRUE)
doMC::registerDoMC(cores = 8)


# ==== METHODS TO COLLECT TISSUE-SPECIFIC EPISPLICED GENES ====
get_p_adj <- function(p_val_list) {
  #' Given a list of a dataframe of p-values, correct the values column-wise for each dataframe
  #'
  #' @param p_val_list a list of dataframe storing the p-values, where columns are samples and rows are genes
  #' @return a list of dataframe storing the FDR-corrected p-values across all genes
  #'

  p_adj_list <- vector("list")
  for (i in seq(length(p_val_list))) {
    gene_id <- p_val_list[[i]]$gene_id
    res_list <- p_val_list[[i]]
    res_list_p <- apply(res_list[, 2:ncol(res_list)], 1, p.adjust, method = "fdr")
    res_list_p <- as.data.frame(cbind(gene_id, t(res_list_p)))
    p_adj_list[[i]] <- res_list_p
  }
  names(p_adj_list) <- histone_type_list
  return(p_adj_list)
}

filter_by_p <- function(r_df_list, p_df_list) {
  #' Given a list of a dataframe of r-values and a list of corresponding p-values, set the unsignificant R-values to 0
  #'
  #' @param r_df_list a list of dataframe storing the r-values, where columns are samples and rows are genes
  #' @param p_df_list a list of dataframe storing the p-values corresponding to those in r_df_list
  #' @return a list of dataframe storing the R-values only for significant correlations
  #'
  
  filtered_df_list <- list()
  for (i in seq(length(r_df_list))) {
    # Select R-/p-values without the gene ID
    r_table <- apply(as.matrix(r_df_list[[i]][, 2:ncol(r_df_list[[i]])]), 2, as.numeric)
    p_table <- apply(as.matrix(p_df_list[[i]][, 2:ncol(p_df_list[[i]])]), 2, as.numeric)

    # Retain R values if p-values < 0.05 (significant correlation)
    r_table[p_table > 0.05] <- NA

    # Complete and add info to R-values dataframe
    r_table <- as.data.frame(r_table)
    r_table <- cbind(r_df_list[[i]][, 1:1], r_table)
    colnames(r_table)[1] <- "gene_id"

    filtered_df_list[[i]] <- r_table
  }
  names(filtered_df_list) <- histone_type_list
  return(filtered_df_list)
}

get_sig_gene_list <- function(r_df_list, method, r_sig = 0.5, p_sig = 0.05, DEU_sig = 2, DHM_sig = 1, 
  compare_threshold = "absolute") {
  #' Get the epispliced genes for each pairwise comparison in different histone context based on either correlation
  #' coefficients, p-values or DEU/DHM significance
  #'
  #' @param r_df_list a list of dataframes storing the R-values (filtered by p-values), where columns are samples and 
  #' rows are genes
  #' @param method "pearcor_r" and "pearcor_p" defines significant genes by R-values or p-values, respectively. "DEU" 
  #' and "his" defines the significant genes only based on DEUs/DHMs events. Use "DEU"/"DHM" for getting the gene sets 
  #' to compare clustering efficiency based on epispliced genes versus DEU/DHM in Subfigure S3.
  #' @param r_sig R values threshold to define epispliced genes (default = 0.5)
  #' @param p_sig p-values threshold to define epispliced genes (default = 0.05)
  #' @param DEU_sig log fold change threshold to define significant DEU genes (default = 2)
  #' @param DHM_sig M-value threshold to define significant DHM genes (default = 1) 
  #' @param compare_threshold "greater"/"lesser" get the values to larger/smaller than the defined threshold 
  #' (one-tailed). The default "absolute" get values with absolute value larger than the threshold (two-tailed)
  #' @return a list of lists of epispliced genes named by pair comparison that they are identified in for each
  #' different histone contexts
  #'

  r_sig_list <- vector("list", length(r_df_list))
  # For each histone mark context
  for (i in seq(length(r_df_list))) {
    r_df <- r_df_list[[i]]

    r_df_sig <- vector("list", ncol(r_df) - 1)
    # For pairwise comparison result array stored in r_df columns (First column is gene names)
    for (j in 2:ncol(r_df)) {
      pairwise_r <- as.numeric(r_df[[j]])

      # Get ids of genes that are significant by comparing the values to correlation R/p-values of DEU/DHM threshold
      if (method == "pearcor_r") {
        if (compare_threshold == "greater") {
          r_df_sig[[j - 1]] <- r_df[pairwise_r >= r_sig & !is.na(pairwise_r), "gene_id"]
        } else if (compare_threshold == "lesser") {
          r_df_sig[[j - 1]] <- r_df[pairwise_r <= r_sig & !is.na(pairwise_r), "gene_id"]
        } else if (compare_threshold == "absolute") {
          r_df_sig[[j - 1]] <- r_df[abs(pairwise_r) >= r_sig & !is.na(pairwise_r), "gene_id"]
        }
      } else if (method == "pearcor_p") {
        r_df_sig[[j - 1]] <- r_df[pairwise_r <= p_sig & !is.na(pairwise_r), "gene_id"]
      } else if (method == "DEU") {
        r_df_sig[[j - 1]] <- r_df[abs(pairwise_r) >= DEU_sig & !is.na(pairwise_r), "gene_id"]
      } else if (method == "DHM") {
        r_df_sig[[j - 1]] <- r_df[abs(pairwise_r) >= DHM_sig & !is.na(pairwise_r), "gene_id"]
      }
    }

    names(r_df_sig) <- colnames(r_df)[2:length(colnames(r_df))]
    r_sig_list[[i]] <- r_df_sig
  }
  return(r_sig_list)
}


get_tissue_specific_indices <- function(r_sig_list) {
  #' Calculate the indices in the list of all pairwise comparisons where a specific tissue would be a part of for all 
  #' tissues
  #'
  #' Since the sets of available tissues in 5 histone contexts are different, the indices, or the columns where a 
  #' specific tissue is found in all pairwise comparisons must be calculated for each context. 
  #' @param r_sig_list a list containing the significant geneset for each histone mark
  #' @return a list containing the set of tissue-specific indices for each histone mark

  all_tissue_index <- vector("list", length(r_sig_list))
  for (i in seq(length(r_sig_list))) {
    all_tissue_name <- names(r_sig_list[[i]])
    tissue_name <- unique(as.vector(sapply(all_tissue_name, function(x) strsplit(x, split = "_")[[1]])))
    tissue_index <- lapply(tissue_name, function(x) grep(x, all_tissue_name))
    names(tissue_index) <- tissue_name
    all_tissue_index[[i]] <- tissue_index
  }
  names(all_tissue_index) <- histone_type_list
  return(all_tissue_index)
}

get_TSEGs <- function(r_sig_list, tissue_spec_ind) {
  #' Pool the Tissue Specific Epispliced Genes (TSEGs) in different histone contexts
  #'
  #' @param r_sig_list a list containing the significant geneset for each histone mark
  #' @param tissue_spec_ind a list containing the set of tissue-specific indices for each histone mark
  #'  

  all_genes_joined <- vector("list")
  for (i in seq(length(r_sig_list))) {
    res <- r_sig_list[[i]]
    tissue_list <- tissue_spec_ind[[i]]
    all_tissues <- lapply(tissue_list, function(x) {
      return(Reduce(union, res[x]))
    })
    names(all_tissues) <- names(tissue_list)
    all_genes_joined[[i]] <- all_tissues
  }
  names(all_genes_joined) <- names(tissue_spec_ind)
  return(all_genes_joined)
}


#==== METHODS TO PREPARE PLOTS & TABLES DATA ====
prep_plot_ecdf <- function(r_df_list) {
  #' Collapse the list of a dataframe containing significant r-values into one dataframe containing the R-values and
  #' histone type where the correlation is found. This method prepares the dataframe used for plot_ecdf()
  #'
  #' @param r_df_list a list of dataframe storing the r-values, where columns are samples and rows are genes
  #' @return a dataframe storing all R-values and the histone mark they associate to. This dataframe is used for 
  #' plotting the ECDF of correlation coefficients in Figure 3A
  #'

  names(r_df_list) <- histone_type_list
  all_r_DHM <- lapply(histone_type_list, function(x) {
    df <- r_df_list[[x]]

    # Extract and store the R-values and histone type
    df <- unlist(as.list(df[, 2:length(df)]))
    df <- as.data.frame(as.numeric(df[!is.na(df)]))
    df$his <- x

    return(df)
  })
  names(all_r_DHM) <- histone_type_list
  all_r_DHM <- do.call(rbind, all_r_DHM)
  colnames(all_r_DHM) <- c("val", "his")
  return(all_r_DHM)
}

plot_ecdf <- function(all_r_compare, dataset_path){
  #' Plot the empirical cumulative distribution for all R-values in Figure 3A
  #' 
  #' @param all_r_compare a dataframe storing all R-values and the histone mark they associate to. This
  #' dataframe could be the concatenate of multiple dataframes generated in prep_plot_ecdf(). The columns contain the
  #' R-values, histone type, and the characteristic distinguishing concatenated dataframes should be named "val", "his"
  #' and "type", respectively. In this analysis, "type"-values were either "p-value > 0.05" or "p-value <= 0.05".
  #' @param dataset_path path to store the resulting plot, preferably path to the specific dataset being analysed
  #' @return none
  #'  
  
  cumulative_r_ecdf <- ggplot(all_r_compare, aes(val, colour = his, linetype = type)) + 
    stat_ecdf(aes(colour = his)) +
    xlab("Pearson correlation R values") + ylab("Probability") +
    scale_linetype_manual(values = c("dotted", "solid"), name = "", 
                          labels = c(expression(p[adj]<=0.05), expression(p[adj]<=1.0))) +
    scale_color_discrete(name = "Histone type") +
    theme_minimal() +
    theme(axis.title = element_text(face = "bold", size = 13),  aspect.ratio = 1,
          panel.border = element_rect(colour = "grey", fill = NA, size = 0.5),
          legend.position = "bottom", legend.box = "vertical", legend.margin = margin()) +
    guides(fill = guide_legend(nrow = 3, byrow = t))

  tiff(filename = paste(dataset_path, "res", paste("cumulative_r_ecdf_", subscript, ".tiff", sep = ""), sep = "/"),
       width = 8, height = 5, units = "in", res = 300)
  print(cumulative_r_ecdf)
  dev.off()
  print(paste("Plotted ECDF for: ", dataset_path, sep = ""))
}

get_TSEGs_no <- function(TSEGs_list) {
  #' Pool the epispliced genes specific for each tissue in different histone contexts (Table 2)
  #'
  #' @param TSEGs_list a list containing Tissue Specific Epispliced Genes (TSEGs) for each histone contexts
  #' @return a list containing the number of TSEGs for each tissue in each histone context
  #'  

  TSEGs_no_df <- lapply(TSEGs_list, function(x) lapply(x, function(y) length(y)))
  TSEGs_no_df <- lapply(TSEGs_no_df, function(x) as.data.frame(cbind(names(x), do.call(rbind, x))))
  TSEGs_no_df <- TSEGs_no_df %>% purrr::reduce(full_join, by = "V1")
  colnames(TSEGs_no_df) <- c("Tissue", histone_type_list)
  TSEGs_no_df$Tissue <- official_name
  return(TSEGs_no_df)
}

get_nonoverlap_TSEG_set <- function(TSEGs_list) {
  #' Pool the epispliced genes specific for each tissue in all histone contexts and discard duplicates 
  #' (Table 2 - last column)
  #'
  #' @param TSEGs_list a list containing Tissue Specific Epispliced Genes (TSEGs) for each histone contexts
  #' @return a list containing the number of TSEGs for each tissue in each histone context
  #'  

  all_DHMtones_TSEGs <- vector("list")

  # For each tissue investigated (the full list of tissues is extracted from names(TSEGs_list[["H3K27ac"]]))
  for (i in seq(length(names(TSEGs_list[["H3K27ac"]])))) {
    tissue <- names(TSEGs_list$H3K27ac)[i]
    all_res_sig <- vector("list")

    # For each histone mark, get the tissue_1-specific epispliced genes from pairwise comparisons tissue_1 involves in
    for (j in seq(length(TSEGs_list))) {
      pairs <- match(tissue, names(TSEGs_list[[j]]))
      all_res_sig[[j]] <- TSEGs_list[[j]][[pairs]]
    }

    # Collapse the sets of TSEGs across all histones for the tissue being considered to 
    gene_set <- Reduce(union, all_res_sig)
    all_DHMtones_TSEGs[[i]] <- gene_set
  }
  return(all_DHMtones_TSEGs)
}

export_no_sig_gene_by_tissue <- function(no_sig_gene_by_tissue, dataset_path) {
  #' Export the number of epispliced genes found for each tissue in each histone context in Table 2
  #' 
  #' @param no_sig_gene_by_tissue a list containing the number of TSEGs for each tissue in each histone context
  #' @param dataset_path path to store the resulting plot, preferably path to the specific dataset being analysed
  #' @return none
  #'  
  
  fwrite(x = no_sig_gene_by_tissue,
        file = paste(dataset_path, "res", paste("no_sig_gene_by_tissue", subscript, ".csv", sep = ""), sep = "/"), 
        quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
}

compute_jaccard <- function(DHM_TSEGs) {
  #' Compute the Jaccard index between the tissues based on the set of tissue_specific TSEGs. The indices are used to
  #' construct the heatmaps in Figure 5 and Subfigure S3
  #'
  #' Jaccard index is calculated for each tissue pair based on the union and overlapping set of TSEGs specific for the
  #' two tissues.
  #' @param DHM_TSEGs a list containing Tissue Specific Epispliced Genes (TSEGs) specific for a histone mark
  #' @return a list containing Jaccard distances between any tissue pair in a specific histone context
  #'  
  
  # Get the intersection between any two TSEGs sets
  distances <- sapply(DHM_TSEGs, function(x) sapply(DHM_TSEGs, function(y) length(intersect(x, y))))

  # Get the union between any two TSEGs sets
  distances_union <- sapply(DHM_TSEGs, function(x) sapply(DHM_TSEGs, function(y) length(union(x, y))))
  diag(distances) <- NA

  # Compute Jaccard index as distance between two tissues
  distances <- distances / distances_union
  return(distances)
}


get_similarity_matrix <- function(TSEGs_list, mat_names) {
  #' Compute the pairwise distance across all tissues based on the sets of TSEGs for each histone context in Figure 5
  #' and Subfigure S3
  #' 
  #' @param TSEGs_list a list containing Tissue Specific Epispliced Genes (TSEGs) for each histone contexts
  #' @param mat_names names list to name the matrices. For the main analysis, the matrices are named by the histone 
  #' types. For the same analysis but with either DEU or DHM genes, only the matrices are named as "DEU"/histone types
  #' @return a list containing Jaccard distances between any tissue pair in each histone context
  #'  
  
  sim_matrices <- vector("list")
  # For each histone mark
  for (i in seq(length(TSEGs_list))) {
    histone_specific_TSEGs <- TSEGs_list[[i]]
    
    # Compute Jaccard index
    distances <- compute_jaccard(DHM_TSEGs = histone_specific_TSEGs) 

    rownames(distances) <- sapply(names(histone_specific_TSEGs), function(x) 
      official_name[grep(x, rownames(epigenomes_annot))[1]])
    colnames(distances) <- rownames(distances)
    
    sim_matrices[[i]] <- distances
  }
  names(sim_matrices) <- mat_names
  print(paste("Highest index: ", max(sapply(sim_matrices, function(x) max(x, na.rm = TRUE)))))
  print(paste("Lowest index: ", min(sapply(sim_matrices, function(x) min(x, na.rm = TRUE)))))

  return(sim_matrices)
}

get_color_scale <- function(sim_matrices, digits = 2) {
  max_distance <- ceiling(max(unlist(lapply(sim_matrices, function(x) max(x, na.rm = TRUE)))) * 10^digits) / 10^digits
  min_distance <- floor(min(unlist(lapply(sim_matrices, function(x) min(x, na.rm = TRUE)))) * 10^digits) / 10^digits
  breaks <- seq(min_distance, max_distance, ((max_distance - min_distance) / 99))
  return(breaks)
}

plot_heatmap <- function(sim_mat, sim_mat_DEU_DHM, dataset_path) {
  #' Plot the heatmap showing the similarity between any tissue pair based on their TSEGs sets in Figure 5 and 
  #' Subfigure S3
  #' 
  #' @param sim_mat a list containing Jaccard distances computed based on TSEGs similarity between any tissue pair in
  #' each histone context. The length of this list is the number of histone marks.
  #' @param sim_mat_DEU_DHM a list containing Jaccard distances computed based on DEU/DHM significant genes-similarity 
  #' between any tissue pair in each histone context. The length of this list is the number of histone marks plus one 
  #' distance matrix for DEU genes
  #' @param dataset_path path to store the resulting plot, preferably path to the specific dataset being analysed
  #' @return none
  #'  
  
  color <- "-RdYlBu2:101"
  epigenomes_colors <- list(
    c("lightgrey", "hotpink", "black"),
    c("plum", "lightblue", "gold", "darkseagreen"),
    c("lightgrey", "hotpink", "black", "MediumPurple"),
    c("plum", "gold")
  )
  
  # Plot similarity heatmaps for tissues based on TSEGs
  breaks <- get_color_scale(sim_matrices = sim_mat)

  tiff(filename = paste(dataset_path, "res", paste("heatmap_", subscript, ".tiff", sep = ""), sep = "/"), 
      width = 16, height = 14, units = "in", res = 200)
  par(mfrow = c(2, 3), mar = c(1, 1, 1, 1) * 10)
  for (i in seq(length(sim_mat))){
    aheatmap(sim_mat[[i]], Rowv = FALSE, Colv = TRUE, scale = "none", main = paste(histone_type_list[[i]]), 
              annCol = epigenomes_annot[names(rownames(sim_mat[[i]])), c("Potency", "Type", "Origin", "Life.stage")],
              annColors = epigenomes_colors, breaks = breaks, color = color, 
              distfun = "maximum", hclustfun = "average", reorderfun = function(d, w) reorder(d, w, agglo.FUN = mean),
              legend = FALSE, annLegend = FALSE, fontsize = 14, cexRow = 0, cexCol = 1, treeheight = 30
              )
  }
  dev.off()

  # Plot similarity heatmaps for tissues based on DEU/DHM genes
  breaks <- get_color_scale(sim_matrices = sim_mat_DEU_DHM)

  tiff(filename = paste(dataset_path, "res", paste("heatmap_", subscript, "_DEU_DHM.tiff", sep = ""), sep = "/"),
      width = 16, height = 14, units = "in", res = 200)
  par(mfrow = c(2, 3), mar = c(1, 1, 1, 1)*10)
  for (i in seq(length(sim_mat_DEU_DHM))){
    aheatmap(sim_mat_DEU_DHM[[i]], Rowv = FALSE, Colv = TRUE, scale = "none", 
              main = paste(c("DEU", paste("DHM", histone_type_list, sep = "-"))[i]), 
              annCol = epigenomes_annot[names(rownames(sim_mat_DEU_DHM[[i]])), c("Potency", "Type", "Origin", "Life.stage")],
              annColors = epigenomes_colors, breaks = breaks, color = color, 
              distfun = "maximum", hclustfun = "average", reorderfun = function(d, w) reorder(d, w, agglo.FUN = mean), 
              legend = FALSE, annLegend = FALSE, fontsize = 14, cexRow = 0, cexCol = 1, treeheight = 30
             )
  }
  dev.off()
}

get_plot_dfs <- function(DEU_df, DHM_dfs, gene, tissue_pair) {
  #' Gather the information used for the correlation computation which will be used by plot_cor() to plot the 
  #' correlation for example genes in Figure 3B
  #' 
  #' Given the dataframes containing either DEU statistics or DHM M-values for all genes, all pairwise comparisons, 
  #' extract the DEU/DHM values between 2 tissues for a specific gene and all histone contexts
  #' @param DEU_df a dataframe of DEU statistics with sample pairs as columns and gene/exons as rows
  #' @param DHM_dfs a list of dataframes containing DHM M-values with sample pairs as columns and 
  #' gene/exons as rows for each histone context
  #' @param gene the gene to be plotted as example
  #' @param tissue_pair the pair of tissues of which DEUs and DHMs are correlated.
  #' @return a list of dataframes containing DEU and DHM values from comparing specific 2 tissues for the example gene
  #' in each histone context
  #'  
  
  return(
    lapply(DHM_dfs, function(x)
      df <- data.frame(
        DEU = DEU_df[[tissue_pair]][DEU_df$gene_id == gene], 
        DHM = x[[tissue_pair]][x$gene_id == gene]) 
    )
  )
}

plot_cor <- function(dfs, plot_title) {
  #' Plot the DEU and DHM values from a comparison for a specific gene and the correlation between them for Figure 3B
  #' 
  #' @param dfs a list of dataframes containing DEU and DHM values from comparing specific 2 tissues for the example 
  #' gene in each histone context. get_plot_dfs() helps creating these dataframes
  #' @param plot_title a string to define the name of the plot. This should contain the gene name and the tissue pairs 
  #' being plotted
  #' @return correlation plots 
  #' 

  name_dfs <- names(dfs)
  plots <- list()
  for (i in seq(length(dfs))){
    df <- dfs[[i]]
    plots[[i]] <- ggscatter(data = df, x = "DEU", y = "DHM", add = "reg.line", 
        add.params = list(color = "#5D2ED2", fill = "lightgray"), conf.int = TRUE,
        font.label = c(8, "plain", "red"), title = paste(plot_title, name_dfs[i], sep = "\n")) +
      stat_cor(label.x = min(df$DEU) + 0.5, label.y = max(df$DHM) + 0.5) +
      ylab(paste("DHM", histone_type_list[[i]], sep = "-")) + xlab("DEU") +
      theme(aspect.ratio = 0.8, title = element_text(size = 9), axis.title = element_text(face = "bold", size = 12))
  }
  names(plots) <- name_dfs
  return(plots)
}

