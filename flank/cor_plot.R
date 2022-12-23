library(ggplot2)
library(reshape2)
library(ggpubr)
library(dplyr)
library(data.table)


# ==== PREPARATION ====
# Define histone list and dataframe
histone_type_list <- c("H3K27ac", "H3K27me3", "H3K36me3", "H3K4me1", "H3K4me3", "H3K9me3")
all_pairs.exp <- as.data.frame(readRDS("all_pairs.exp.RDS"))
all_pairs.his <- readRDS("all_pairs.his_list.RDS")
all_pairs.his <- lapply(all_pairs.his, as.data.frame)

pearson_r <- readRDS("all_res_list.pearcor_r.RDS") # TODO: ????
pearson_p_adj <- readRDS("all_res_list.pearcor_padj.RDS")


# ==== FIGURE 3A: PLOT CORRELATION CUMULATIVE DISTRIBUTION ====
filter_by_p <- function(r_df, p_df) { 
  #' Filter R-values by FDR-adjusted p-values
  #' 
  #' From the list of R-values for all comparisons, keep only those with FDR-adjusted p-values < 0.05
  #' @param r_df list of dataframe containg R-values, where columns are tissue pair and rows are genes
  #' @param p_df list of dataframe containg FDR-adjusted p-values, where columns are tissue pair and rows are genes
  #' 
  r_df_filtered <- list()

  for (i in 1:length(histone_type_list)) {
    r_table <- as.matrix(r_df[[i]][, 2:ncol(r_df[[i]])])
    p_table <- as.matrix(p_df[[i]][, 2:ncol(p_df[[i]])])
    r_table[p_table > 0.05] <- NA
    r_table <- as.data.frame(r_table)
    r_table <- cbind(r_df[[i]][, 1:2], r_table)
    r_df_filtered[[i]] <- r_table
  }

  return(r_df_filtered)
}

# Filter for p-values with FDR < 0.05
pearson_r_filtered <- filter_by_p(pearson_r, pearson_p_adj)


get_occurence <- function(r_df, significance) {
  #' Prepare a table of occurence for R-values
  #' 
  #' @param r_df list of dataframe containg R-values, where columns are tissue pair and rows are genes
  #' @param significance the significance level used to filter the R-values in r_df
  #' 
  all_r_vals <- list()
  for (i in 1:length(histone_type_list)) {
    occ_df <- melt(r_df[[i]][, 2:ncol(r_df[[i]])])
    occ_df$H_type <- histone_type_list[i]
    all_r_vals[[i]] <- occ_df
  }
  all_r_vals <- do.call(rbind, all_r_vals)
  all_r_vals <- all_r_vals[!is.na(all_r_vals$value), ]
  all_r_vals$significance <- significance
  return(all_r_vals)
}

# Get occurence of R-values for unfiltered result and for R-values with FDR < 0.05
pearson_r_occ_filtered <- get_occurence(pearson_r_filtered, significance = 0.05)
pearson_r_occ_unfiltered <- get_occurence(pearson_r, significance = 1.0)
all_pearson_r_occ <- rbind(pearson_r_occ_filtered, pearson_r_occ_unfiltered)

# Plot ECDF for R-values
cumulative_r_ecdf <- ggplot(all_pearson_r_occ, aes(value, colour = H_type, linetype = significance)) +
  stat_ecdf(aes(colour = H_type)) +
  xlab("Pearson correlation R values") +
  ylab("Probability") +
  scale_linetype_manual(values = c("solid", "dotted")) +
  theme_bw() +
  theme(aspect.ratio = 1)


# ==== FIGURE 3B: PLOT CORRELATION FOR EXAMPLE GENES ====
prepare_cor_array <- function(gene, tissue_pair, H_index) {
  #' Prepare dataframe for correlation plot
  #'
  #' Return a dataframe of 2 columns for Exon usage and M-value of a gene
  #' @param gene The gene with Exon usage and M-value to be plotted
  #' @param tissue_pair Combination of tissues, separated by '_'
  #' @param H_index Index of histone type
  #'
  exp <- all_pairs.exp[all_pairs.exp$gene_id == gene, tissue_pair]
  tissue_pair <- gsub("trophoblastcell", "trophoblast", tissue_pair)
  his <- all_pairs.his[[H_index]]
  his <- his[his$gene_id == gene, tissue_pair]
  exp_his <- as.data.frame(cbind(exp, his))
  return(exp_his)
}

make_cor_plot <- function(df, stat = TRUE) {
  #' Plot DHMs and DEUs and the correlation between them
  #'
  #' The DEU and DHM values of a gene specific for a tissue comparison and a histone type are plotted with correlation
  #' @param df dataframe of 2 columns DEU and DHM from prepare_cor_array()
  #' @param stat include correlation plot or not
  #' 
  if (stat == TRUE) {
    plot <- ggplot(df, aes(x = exp, y = his)) +
      geom_point() +
      geom_smooth(method = "lm", formula = y ~ x) +
      stat_cor(method = "pearson", label.x = 0, label.y = -1.5, digits = 2)
  } else {
    plot <- ggplot(df, aes(x = exp, y = his)) +
      geom_point() +
      geom_smooth(method = "lm", formula = y ~ x) 
      }
  plot <- plot +
    theme_bw() +
    theme(
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      panel.background = element_blank(), aspect.ratio = 1
    ) +
    xlab(NULL) + ylab(NULL)
  return(plot)
}

make_cor_plot_example <- function(gene, tissue_pair, skipping_cor) {
  #' Plot DHMs and DEUs and the correlation between them
  #'
  #' Generic function to plot correlation. The gene and tissue pair to be plotted are defined here for make_cor_plot()
  #' @param gene gene of interest
  #' @param tissue_pair tissue pair to compare DEU and DHM 
  #' @param skipping_cor list of index of histone type to be skipped
  #' 
  plots <- list()
  for (i in 1:length(histone_type_list)) {
    cor_array <- prepare_cor_array(gene, tissue_pair, i)
    if (ncol(cor_array) == 2) {
      colnames(cor_array) <- c("exp", "his")

      if (i %in% skipping_cor) {
        plot <- make_cor_plot(cor_array, stat = FALSE)
      } else {
        plot <- make_cor_plot(cor_array, stat = TRUE)
      }
      plots[[i]] <- plot
    } else {
      print(i)
      plots[[i]] <- NULL
    }
  }
  return(plots)
}
FGFR2_plot <- make_cor_plot_example("LMNB1", "neuronalstemcell_pancreas", seq(2, 6))
LMNB1_plot <- make_cor_plot_example("FGFR2", "mesenchymalstemcell_sigmoidcolon", c(2, 4, 6))


# ==== COMBINE ALL PLOTS TOGETHER ====
example_gene_plots <- ggarrange(
  LMNB1_plot[[1]], LMNB1_plot[[3]], LMNB1_plot[[5]], FGFR2_plot[[2]],
  ncol = 2, nrow = 2, align = "hv"
)
example_gene_plots <- annotate_figure(example_gene_plots,
  bottom = text_grob("Differential exon usage (Exon dispersion statistic)"),
  left = text_grob("Differential histone modification (M-value)", rot = 90)
)
CorPlot <- ggarrange(example_gene_plots, cumulative_r_ecdf,
  nrow = 1, ncol = 2, align = "hv"
)
tiff("figure3.tiff", width = 15, height = 12, units = "in", res = 200) # save pdf 20*8
CorPlot
dev.off()
