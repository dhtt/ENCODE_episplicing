## ---------------------------
##
## Script name: analysis_pipeline.R
##
## Purpose of script: Perform the main analysis on DEU-DHM correlation data 
##
## Author: Trang Do
##
## Copyright (c) Trang Do, 2018
## Email: dhttrang@bioinformatik.uni-saarland.de
##
## ---------------------------
##
## Notes:
## - Description: This script calls the analysis methods in analysis_methods.R to extract the epispliced genes to 
##   fulfill the main aim of the study. Figures and tables along the manuscript are generated from this script, 
##   including Figure 3, Figure 5, Figure 6, Subfigure S3, Subfigure S4, Subfigure S1B, Table 2
## - Preceeding script: get_cor.R
## - Succeeding script: -
##
## ---------------------------

# ==== LOAD PACKAGES =====
library("reshape2", quietly = TRUE)
library("data.table", quietly = TRUE)
library("dplyr", quietly = TRUE)
library("questionr", quietly = TRUE)
library("ggrepel", quietly = TRUE)
library("scales", quietly = TRUE)
library("org.Hs.eg.db", quietly = TRUE)
library("clusterProfiler", quietly = TRUE)
library("ggforce", quietly = TRUE)
library("GO.db", quietly = TRUE)
library("tidyverse", quietly = TRUE)
library("tidytext", quietly = TRUE)
library("ggpubr", quietly = TRUE)
library("optparse", quietly = TRUE)

# ==== PREPARATION ====
# Parse arguments for input/output settings
option_list <- list(
  make_option("--dataset_path",
    type = "character",
    help = "path to the folder where the general results from the analysis are stored for a specific dataset",
    metavar = "character",
    default = "general_analysis_results/datasets/True_Pearson"
  )
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)


histone_type_list <- c("H3K27ac", "H3K27me3", "H3K36me3", "H3K4me3", "H3K9me3")
DEU_df_path <- paste(opt$dataset_path, "DEU_df.RDS", sep = "/")
DHM_dfs_path <- paste(opt$dataset_path, "DEU_df.RDS", sep = "/")
genewise_odds_ratios_plot_path <- paste(opt$dataset_path, "res", "genewise_or.tiff", sep = "/")
odds_ratios_plot_path <- paste(opt$dataset_path, "res", "groupwise_or.tiff", sep = "/")

# ==== COMPUTE GLOBAL ODDS RATIO ====
# ---- Methods ----
custom_reduce <- function(df_col){
  return(Reduce(function(a, b) a|b, as.logical(df_col)))
}

prepare_global_overlap_OR <- function(DEU_df, DHM_dfs){
  #' Prepare a dataframe with DHM occurence and DEU occurrence from a retrieved from get_cor()
  #'
  #' This function reads the DHMs results across all pairwise comparisons and summarizes the DHM events for each exon
  #' by marking the DHM = TRUE if it occurs in at least 1 comparison. DEU results are collapsed in this manner, too.
  #' @param DEU_df a dataframe containing DEU statistics filtered by first exons and significant DEUs
  #' @param DHM_dfs a list of dataframes containing DHM M-values filtered by first exons and significant DHMs for 
  #' each histone context
  #' @return a dataframe with DEU and DHM events marked for each exon
  #'
  
  # For each histone context and each exon, if DHM occurs in any pairwise comparison, mark for that exon DHM = TRUE
  DHM_any_pair <- lapply(DHM_dfs, function(x) {
    # Convert M-values into boolean values, where DHM = TRUE if absolute M-values > 0 and FALSE otherwise
    df = abs(x[, 3:ncol(x)]) > 0
    df_ = apply(df, 1, custom_reduce)
  })
  DHM_any_pair_df <- do.call("cbind", DHM_any_pair)

  # For each exon, if DHM occurs in any histone context, mark for that exon DHM = TRUE
  DHM_any_histone <- as.data.frame(apply(DHM_any_pair_df, 1, custom_reduce))

  # Prepare DEU data
  DEU_any_pair <- DEU_df[, 3:ncol(DEU_df)]
  # Convert DEU into boolean values, where DEU = TRUE if DEU > 0 and FALSE otherwise
  DEU_any_pair <- DEU_any_pair > 0
  # For each exon, if DEU occurs in any pairwise comparison, mark for that exon DEU = TRUE
  DEU_any_pair <- as.data.frame(apply(DEU_any_pair, 1, custom_reduce))

  # Prepare a dataframe storing both DHMs and DEU
  DHM_DEU_df_boolean <- do.call('cbind', list(DEU_df$gene_id, DHM_any_pair_df, DHM_any_histone, DEU_any_pair))

  colnames(DHM_DEU_df_boolean) = c("gene_id", histone_type_list, "dhm", "deu")
  return(DHM_DEU_df_boolean)
}

get_global_overlap_OR <- function(DHM_DEU_df_boolean){
  #' Compute the Odds Ratio between any type of DHM and DEU to measure the global co-occurrence (NOT gene-specific).
  #' The results from this function are shown in Table 1
  #'
  #' @param DHM_DEU_df_boolean a dataframe with DEU and DHM events marked for each exon
  #' @return a list of DHM-DEU contigency tables for each histone context
  #'
  
  DEU <- DHM_DEU_df_boolean[[ncol(DHM_DEU_df_boolean)]]

  DEU_DHM_table_list <- list()
  for (i in seq(2, ncol(DHM_DEU_df_boolean) - 1)){
    DHM <- DHM_DEU_df_boolean[[i]]
    DEU_DHM_table <- table(DHM, DEU)
    print(DEU_DHM_table)
    print(odds.ratio(DEU_DHM_table))
    DEU_DHM_table_list[[i]] <- DEU_DHM_table
  }
  return(DEU_DHM_table_list)
}

# ---- Compute global odds ratios (Table 1) ----
DEU_df <- readRDS(DEU_df_path)
DHM_dfs <- readRDS(DHM_dfs_path)
colnames(DEU_df) = gsub("trophoblastcell", "trophoblast", colnames(DEU_df))

DHM_DEU_df_boolean <- prepare_global_overlap_OR(DEU_df = DEU_df, DHM_dfs = DHM_dfs)
DEU_DHM_table_list <- get_global_overlap_OR(DHM_DEU_df_boolean)

 
# ==== COMPUTE ODDS RATIO FOR EACH GENES ====
# ---- Methods ----
get_odds_ratio_each_genes <- function(DHM_DEU_df_boolean){
  #' Compute the Odds Ratio between any type of DHM and DEU for a specific gene.
  #'
  #' @param DHM_DEU_df_boolean a dataframe with DEU and DHM events marked for each exon
  #' @return a list containing a dataframe with genewise odd ratios and a dataframe with associated p-values for 
  #' all histone context
  #'
  
  # ---- Compute genewise odds ratio ----
  or_df <- DHM_DEU_df_boolean %>%
    dplyr::group_by(gene_id) %>%
    dplyr::summarise(H3K27ac = ifelse(length(unique(H3K27ac)) == 2 && length(unique(deu)) == 2,
                                      yes = odds.ratio(as.factor(H3K27ac), as.factor(deu))$OR, no = NA),
                     H3K27me3 = ifelse(length(unique(H3K27me3)) == 2 && length(unique(deu)) == 2,
                                       yes = odds.ratio(as.factor(H3K27me3), as.factor(deu))$OR, no = NA),
                     H3K36me3 = ifelse(length(unique(H3K36me3)) == 2 && length(unique(deu)) == 2,
                                       yes = odds.ratio(as.factor(H3K36me3), as.factor(deu))$OR, no = NA),
                     H3K4me3 = ifelse(length(unique(H3K4me3)) == 2 && length(unique(deu)) == 2,
                                      yes = odds.ratio(as.factor(H3K4me3), as.factor(deu))$OR, no = NA),
                     H3K9me3 = ifelse(length(unique(H3K9me3)) == 2 && length(unique(deu)) == 2,
                                      yes = odds.ratio(as.factor(H3K9me3), as.factor(deu))$OR, no = NA),
                     dhm = ifelse(length(unique(dhm)) == 2 && length(unique(deu)) == 2,
                                  yes = odds.ratio(as.factor(dhm), as.factor(deu))$OR, no = NA),
                     n_exons = n()
                     ) %>%
    as.data.frame()
  rownames(or_df) = or_df$gene_id 
  or_df <- dplyr::select(or_df, -gene_id)
  
  
  # ---- Compute genewise odds ratio p-values ----
  p_df <- DHM_DEU_df_boolean %>%
    dplyr::group_by(gene_id) %>%
    dplyr::summarise(H3K27ac = ifelse(length(unique(H3K27ac)) == 2 && length(unique(deu)) == 2,
                                      yes = odds.ratio(as.factor(H3K27ac), as.factor(deu))$p, no = NA),
                     H3K27me3 = ifelse(length(unique(H3K27me3)) == 2 && length(unique(deu)) == 2,
                                       yes = odds.ratio(as.factor(H3K27me3), as.factor(deu))$p, no = NA),
                     H3K36me3 = ifelse(length(unique(H3K36me3)) == 2 && length(unique(deu)) == 2,
                                       yes = odds.ratio(as.factor(H3K36me3), as.factor(deu))$p, no = NA),
                     H3K4me3 = ifelse(length(unique(H3K4me3)) == 2 && length(unique(deu)) == 2,
                                      yes = odds.ratio(as.factor(H3K4me3), as.factor(deu))$p, no = NA),
                     H3K9me3 = ifelse(length(unique(H3K9me3)) == 2 && length(unique(deu)) == 2,
                                      yes = odds.ratio(as.factor(H3K9me3), as.factor(deu))$p, no = NA) ,
                     dhm = ifelse(length(unique(dhm)) == 2 && length(unique(deu)) == 2,
                                  yes = odds.ratio(as.factor(dhm), as.factor(deu))$p, no = NA),
                     n_exons = n()
                     ) %>%
    as.data.frame()
  rownames(p_df) = p_df$gene_id 
  p_df <- dplyr::select(p_df, -gene_id)

  return(list(or_df, p_df))
}

prepare_genewise_OR <- function(odds_ratio_result){
  #' Prepare genewise odds ratios for plotting.
  #'
  #' @param odds_ratio_result a list containing a dataframe with genewise odd ratios and a dataframe with associated 
  #' p-values for all histone context
  #' @return a cleaned up dataframe for plotting genewise odds ratios containing odds ratio, p-value, number of exons
  #' for each gene in each histone context
  
  rownames(odds_ratio_result[[2]]) = rownames(odds_ratio_result[[1]])

  odds_ratio_result_melt <- lapply(odds_ratio_result, function(x) {
    df <- as.data.frame(x)
    colnames(df) = c(histone_type_list, "dhm", "n_exons")
    df$gene = rownames(x)
    df = reshape2::melt(df, id.vars=c('gene', 'n_exons'), value.name = "val", variable.name = "his", na.rm = F)
    return(df)})

  odds_ratio_melt_df <- cbind(odds_ratio_result_melt[[1]], odds_ratio_result_melt[[2]]$val)
  odds_ratio_melt_df$his = gsub("dhm", "All histones", odds_ratio_melt_df$his)
  colnames(odds_ratio_melt_df) = c("gene_id", "n_exons", "his", "val", "p_val")
  odds_ratio_melt_df = odds_ratio_melt_df %>%
    dplyr::group_by(his, n_exons) %>%
    mutate(
      p_adj = p.adjust(p_val, "fdr"),
      sig = ifelse(p_adj <= 0.05, TRUE, FALSE)
    ) %>%
    ungroup()

  return(odds_ratio_melt_df)
}

squish_trans <- function(from, to, factor) {
  #' Method adopted from sysilviakim/Kmisc to scale x axis
  #' 
  trans <- function(x) {
    
    if (any(is.na(x))) return(x)
    
    # get indices for the relevant regions
    isq <- x > from & x < to
    ito <- x >= to
    
    # apply transformation
    x[isq] <- from + (x[isq] - from)/factor
    x[ito] <- from + (to - from)/factor + (x[ito] - to)
    
    return(x)
  }
  
  inv <- function(x) {
    
    if (any(is.na(x))) return(x)
    
    # get indices for the relevant regions
    isq <- x > from & x < from + (to - from)/factor
    ito <- x >= from + (to - from)/factor
    
    # apply transformation
    x[isq] <- from + (x[isq] - from) * factor
    x[ito] <- to + (x[ito] - (from + (to - from)/factor))
    
    return(x)
  }
  
  # return the transformation
  return(trans_new("squished", trans, inv))
}

plot_genewise_ORs <- function(odds_ratio_melt_df, filename){
  #' Plot the genewise odds ratios for all genes in Subfigure S2
  #'
  #' @param odds_ratio_melt_df a cleaned up dataframe for plotting genewise odds ratios with odds ratio, p-value, 
  #' number of exons for each gene in each histone context
  #' @param filename path to file where the plot should be saved
  #' @return genewise odds ratio plot
  #'
  
  ORs_plot <- ggplot(odds_ratio_melt_df[is.finite(odds_ratio_melt_df$val), ], 
        aes(x = as.numeric(as.character(n_exons)), y = val, group = n_exons, label = gene_id)) +
    geom_col(position = "dodge2", col = "#472d7b", fill = "#472d7b") +
    geom_point(size = 3, show.legend = T, aes(col = sig),
              alpha = if_else(odds_ratio_melt_df$sig[is.finite(odds_ratio_melt_df$val)] == TRUE &  
                              odds_ratio_melt_df$val[is.finite(odds_ratio_melt_df$val)] > 1, 0.5, 0)) +
    geom_label_repel(aes(label = ifelse(sig == TRUE & val > 1, as.character(gene_id), "")), box.padding = 1, 
                    max.overlaps = Inf, size = 5, max.time = 2, segment.color = "grey50", nudge_x = 30, nudge_y = 10) +
    theme_minimal() +
    facet_wrap(facets = vars(his), nrow = 3) +
    ylab("Odd ratio") + xlab("Number of exons per gene") + 
    scale_colour_manual(labels = c("Not significant", expression(paste(OR>1, " & ", p[FET]<=0.05))), 
                        values = c("#5D2ED2", "red"), name = "") +
    scale_y_continuous(breaks = seq(0, 60, 10)) +
    scale_x_continuous(trans = squish_trans(125, 190, 3), breaks = c(seq(0, 125, 25), seq(160, 200, 30))) +
    theme(legend.position = "top", legend.text = element_text(size = 15), aspect.ratio = 0.6,
          axis.title = element_text(size = 15), axis.title.x = element_text(vjust = 0), 
          axis.text = element_text(size = 12), strip.background = element_rect(fill = "#5D2ED2"),
          panel.border = element_rect(fill = NA, size = 0.25), panel.grid.minor.y = element_blank(),
          strip.text = element_text(face = "bold", size = 14, colour = "white")
          ) 

  tiff(filename, width = 12, height = 13, res = 300, units = "in")
  print(ORs_plot)
  dev.off()

  return(ORs_plot)
}

plot_genewise_ORs_boxplot <- function(odds_ratio_melt_df){
  #' Plot the genewise odds ratios and their distribution in each histone context for all genes in Figure 2A
  #'
  #' @param odds_ratio_melt_df a cleaned up dataframe for plotting genewise odds ratios with odds ratio, p-value, 
  #' number of exons for each gene in each histone context
  #' @return genewise odds ratio plot
  #'
  
  cumulative_odd_ecdf_box <- ggplot(odds_ratio_melt_df, aes(val)) +   
    geom_hline(yintercept = 1, color= "red", size = 0.4) + 
    geom_boxplot(aes(x = his, y = val), outlier.size = 0.8, col = "#5D2ED2", fill = "lightgrey", alpha = 0.1) +
    geom_point(
      data = odds_ratio_melt_df[odds_ratio_melt_df$sig == TRUE & odds_ratio_melt_df$val > 1 & is.finite(odds_ratio_melt_df$val), ], 
      aes(y = val, x = his), col = "red", alpha = 0.5, size = 2) +
    facet_zoom(ylim = c(0, 10), zoom.size = 1, show.area = F)+ 
    ylab("Odd ratio") + xlab("") +
    scale_y_continuous(breaks = c(0, 1, seq(10, 70, 10))) +
    scale_color_viridis_d(end = 0.9, name = "", direction = 1) +
    scale_fill_viridis_d(end = 0.9, name = "", direction = 1) +
    theme_light() + 
    theme(axis.title.x = element_text(face = "bold"), axis.title.y = element_text(face = "bold", size = 14),
          axis.text.x = element_text(hjust = 1, angle = 30, size = 13),
          axis.text.y = element_text(vjust = c(1, 0, rep(1, 6)), face = c("plain", "bold", rep("plain", 6)),
                                     colour = c("black", "red", rep("black", 6)), size = 12),
          axis.ticks.y = element_line(colour = c("black", "red", rep("black", 6)), size = c(0.5, 0.4, rep(0.5, 6))),
          legend.position = "bottom", panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = "grey", fill = NA, size = 0.5))

  return(cumulative_odd_ecdf_box)
}

# ---- Prepare data for plotting genewise ORs (Subfigure S2 & Figure 2A) ----
odds_ratio_res <- get_odds_ratio_each_genes(DHM_DEU_df_boolean)
odds_ratio_melt_df <- prepare_genewise_OR(odds_ratio_result = odds_ratio_res)

# ---- Plot the genewise ORs (Subfigure S2 & Figure 2A) ----
plot_genewise_ORs(odds_ratio_melt_df = odds_ratio_melt_df, filename = genewise_odds_ratios_plot_path)
cumulative_odd_ecdf_box <- plot_genewise_ORs_boxplot(odds_ratio_melt_df = odds_ratio_melt_df)

# ==== COMPUTE ODDS RATIO FOR GROUPS ====
# ---- Methods ----
get_entrez_id <- function(gene_list){
  #' Given a list of HGNC gene IDs, return the GO terms annotation associated to these genes 
  #'
  #' @param gene_list a list of HGNC gene ids
  #' @return a dataframe containing the GO terms annotations of given genes
  #'
  
  entrez_gene_list <- AnnotationDbi::select(org.Hs.eg.db, keys = gene_list, columns = c("ENTREZID", "SYMBOL"), 
                                            keytype = "SYMBOL")
  entrez_gene_list <- entrez_gene_list[[2]]
  entrez_gene_list <- enrichGO(entrez_gene_list, OrgDb="org.Hs.eg.db", ont = "BP", qvalueCutoff = 0.01, 
                               pAdjustMethod = "fdr", readable = TRUE)
  entrez_gene_list <- entrez_gene_list@result

  return(entrez_gene_list)
}

get_GO_terms_level3 <- function(){
  #'  Prepare a dataframe storing all the level 2 and level 3 terms in the GO terms hierachy 
  #'
  #' @return a dataframe with the GO terms at level 2 and level 3 in the GO terms hierachy for all genes possible and
  #' the genes associated to those universal terms
  #'
  
  # Get all terms related to biological processes involved in development. These are children of the term GO"0008150 
  # excluding "GO:0006791", "GO:0006794", "GO:0015976"
  BP_terms <- AnnotationDbi::select(GO.db, columns = c("TERM"), keytype="GOID",
            keys = setdiff(unlist(as.list(GOBPCHILDREN["GO:0008150"])), c("GO:0006791", "GO:0006794", "GO:0015976")))
            
  # Get terms at level 3 GO terms hierarchy. Those are children of the BP_terms above 
  BP_terms <- lapply(BP_terms$GOID, function(x) data.frame(level_2 = x, level_lv3 = unlist(as.list(GOBPCHILDREN[x]))))
  BP_terms <- do.call(rbind, BP_terms)
  BP_terms <- BP_terms[!is.na(BP_terms$level_lv3), ]

  # Add annotation to the level 2 BP_terms 
  BP_terms$level_2_annot <- AnnotationDbi::select(GO.db, keys = BP_terms$level_2, columns = "TERM")$TERM

  # Create dataframe with only level 2 annots, level 3 annots, and level 3 GO terms 
  BP_terms_lv3_df <- cbind(AnnotationDbi::select(GO.db, keys = BP_terms$level_lv3, columns = c("TERM", "GOID")), 
                          BP_terms$level_2_annot)
  colnames(BP_terms_lv3_df) <- c("GOID", "TERM", "level2_TERM")
  BP_terms_lv3_df <- BP_terms_lv3_df %>% dplyr::filter(!is.na(TERM))

  return(BP_terms_lv3_df)
}

get_enrichment_ratio_level3 <- function(annotated_terms_df, all_GO_term_df, DHM_DEU_df_boolean, his_type){
  #' Compute the odds ratio as enrichment score for each of the level 3 GO terms from the list of epispliced genes 
  #' associated to this term, specifically for a histone context
  #'
  #' @param annotated_terms_df a dataframe containing the GO terms annotations of given genes
  #' @param all_GO_term_df a dataframe with the GO terms at level 2 and level 3 in the GO terms hierachy for all genes
  #' possible and the genes associated to those universal terms
  #' @param DHM_DEU_df_boolean a dataframe with DEU and DHM events marked for each exon
  #' @param his_type histone type being assessed
  #' @return a dataframe for enrichment score for each of the level 3 GO terms based on odds ratios of epispliced genes
  #' for the histone being assessed
  #'
  
  enrich_terms <- lapply(all_GO_term_df$GOID, function(x){
    # Get the genes associated to the specific term
    gene_list <- unique(unlist(lapply(
      annotated_terms_df$geneID[annotated_terms_df$ID %in% unlist(as.list(GOBPOFFSPRING[[x]]))], 
      function(y) strsplit(y, "/")[[1]])))
    
    # Look up the DEU/DHM occurrences for the term-associated genes
    cont_df <- DHM_DEU_df_boolean[DHM_DEU_df_boolean$gene_id %in% gene_list, ] %>% dplyr::select(his_type, deu)
    
    # Compute odds ratios and p-values from the DEU-DHM co-occurrence in this group of genes
    if (dim(cont_df)[1] == 0) {
      # Skip computing if no genes for the term overlap with epispliced genes
      cont_table <- c(NA, NA, length(gene_list))
    }
    else {
      cont_table <- table(cont_df)
      if (ncol(cont_table) == 1 | nrow(cont_table) == 1){
        # Skip computing if contingency table is not available
        cont_table <- c(NA, NA, length(gene_list))
      }
      else {
        cont_table <- apply(cont_table, 2, as.numeric)
        odd <- odds.ratio(cont_table)
        cont_table <- c(round(odd$OR, 3), odd$p, length(gene_list))
      }
    }
    return(cont_table)
  })

  # Convert the list of results into a dataframe
  enrich_terms <- as.data.frame(t(as.data.frame(enrich_terms)))
  colnames(enrich_terms) <- c("enrichment", "p_val", "no_genes")
  enrich_terms <- enrich_terms %>%
    dplyr::mutate(term = all_GO_term_df$level2_TERM,
                  term_spec = all_GO_term_df$TERM)
  return(enrich_terms)
}

get_enrichment_ratio_level3_his <- function(annotated_terms_df, all_GO_term_df, DHM_DEU_df_boolean){
  #' Compute the odds ratio as enrichment score for each of the level 3 GO terms from the list of epispliced genes 
  #' associated to this term for each histone context
  #' 
  #' @param annotated_terms_df a dataframe containing the GO terms annotations of given genes
  #' @param all_GO_term_df a dataframe with the GO terms at level 2 and level 3 in the GO terms hierachy for all genes
  #' possible and the genes associated to those universal terms
  #' @param DHM_DEU_df_boolean a dataframe with DEU and DHM events marked for each exon
  #' @return a dataframe for enrichment score for each of the level 3 GO terms based on odds ratios of epispliced genes
  #' for all histone 


  all_res <- list()
  for (i in seq(length(c(histone_type_list, "dhm")))){
    his <- c(histone_type_list, "dhm")[i]
    print(his) 

    # Compute the odds ratio as enrichment score for each of the level 3 GO terms for a histone mark
    enrich_df <- get_enrichment_ratio_level3(annotated_terms_df = annotated_terms_df, all_GO_term_df = all_GO_term_df, 
      DHM_DEU_df_boolean = DHM_DEU_df_boolean, his_type = his)
    
    # Compute the baseline odds ratio from all epispliced genes (not from a term-specific group of genes)
    enrich_baseline <- table(DHM_DEU_df_boolean[[his]], DHM_DEU_df_boolean[["deu"]])
    enrich_baseline <- apply(enrich_baseline, 1, as.numeric)
    enrich_baseline_OR <- odds.ratio(enrich_baseline)
    enrich_baseline <- data.frame(enrichment = round(enrich_baseline_OR$OR, 3), p_val = enrich_baseline_OR$p, 
                                  no_genes = length(unique(DHM_DEU_df_boolean$gene_id)),
                                  term = "BASELINE", term_spec="BASELINE", row.names = "BASELINE")

    # Combine the odds ratios for the gene groups annotated to a specific term and the odds ratios for all episplced 
    # genes                              
    enrich_df <- rbind(enrich_df, enrich_baseline)
    enrich_df$p_adj <- p.adjust(enrich_df$p_val, "fdr")
    enrich_df$his <- his
    enrich_df$compare_baseline <- enrich_df$enrichment > enrich_baseline$enrichment
    all_res[[i]] <- enrich_df
  }
  names(all_res) <- c(histone_type_list, "dhm")
  return(all_res)
}

plot_groupwise_ORs <- function(enrich_his_lv3_df){
  #' Plot the odds ratios for the gene groups annotated to a specific level 3 GO biological term against the odds 
  #' ratios for all epispliced genes to identify which terms are more enriched in DEU-DHM co-occurrence context.
  #' The plot is shown as Figure 2B
  #'
  #' @param enrich_his_lv3_df  dataframe for enrichment score for each of the level 3 GO terms based on odds ratios of
  #' epispliced genes for all histone 
  #' @return a plot for groupwise odds ratios 
  #'

  odd_ratio_annot <- ggplot(enrich_his_lv3_df, 
                         aes(y = enrichment, x = reorder_within(term_spec, -enrichment, his), col = term, fill = term)) +
  scale_x_reordered() + facet_wrap(. ~ his, scales = "free_x", ncol = 3) +
  geom_segment(aes(x = reorder_within(term_spec, -enrichment, his), xend = reorder_within(term_spec, -enrichment, his),
                   y = 0, yend = enrichment),
               size = 0.75, color = ifelse(enrich_his_lv3_df$term == "BASELINE", "red", "grey")) + 
  geom_point(aes(shape=compare_baseline), size = 4) +
  geom_point(aes(shape=compare_baseline), col = "white", alpha = 0.75, size = 1) +
  theme_minimal() +
  ylab("Odds ratio") + ylim(c(0, 8.5)) + 
  scale_fill_discrete(name = "Level 2 GO terms for Biological Functions") + 
  scale_color_discrete(name = "Level 2 GO terms for Biological Functions") + 
  scale_shape_manual(values=c(18, 19), labels = c("< BASELINE", "> BASELINE"), name="")+
  theme(legend.position="bottom", legend.text = element_text(size = 13), 
        legend.title = element_text(size = 14, face = "bold"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 12, colour = "black"),
        axis.title.x = element_blank(), axis.title.y = element_text(face = "bold", size = 14),
        axis.text.y = element_text(size = 13),
        strip.text = element_text(face = "bold", size = 12, colour = "white"),
        strip.background = element_rect(colour = "#5D2ED2", fill = "#5D2ED2"),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "#5D2ED2", fill = NA, size = 0.25)) +
  guides(fill = guide_legend(title.position = "top", title.hjust = 0.5, ncol = 3, byrow = F),
         shape = guide_legend(title.position = "top", title.hjust = 0.5, ncol = 1, byrow = F)) 

  return(odd_ratio_annot)
}

# ---- Compute the groupwise ORs (Figure 2B) ----
all_GO_annotation <- get_entrez_id(gene_list = as.character(unique(odds_ratio_melt_df$gene_id)))
BP_terms_lv3_df <- get_GO_terms_level3()
enrich_his_lv3 <- get_enrichment_ratio_level3_his(annotated_terms_df = all_GO_annotation, all_GO_term_df = BP_terms_lv3_df, 
                                                  DHM_DEU_df_boolean = DHM_DEU_df_boolean)
enrich_his_lv3_df <- do.call(rbind, enrich_his_lv3) %>% 
  dplyr::filter(!is.na(p_adj) & p_adj < 0.05) %>%
  dplyr::group_by(his) %>%
  dplyr::top_n(15, -p_adj) %>%
  dplyr::top_n(15, enrichment) %>%
  dplyr::mutate(term = str_to_title(term), term_spec = str_to_title(term_spec), his = gsub("dhm", "All histones", his),
                term = gsub("Baseline", "BASELINE", term), term_spec = gsub("Baseline", "BASELINE", term_spec))
enrich_his_lv3_df$term_spec <- stringr::str_to_sentence(enrich_his_lv3_df$term_spec)
enrich_his_lv3_df$term_spec <- gsub("Baseline", "BASELINE", enrich_his_lv3_df$term_spec)
odd_ratio_annot <- plot_groupwise_ORs(enrich_his_lv3_df)

# ---- Plot the genewise and groupwise ORs (Figure 2A & 2B) ----
tiff(odds_ratios_plot_path, res = 300, units = "in", width = 12, height = 14)
ggarrange(plotlist = list(cumulative_odd_ecdf_box, odd_ratio_annot), nrow = 2, labels = c("A", "B"), 
          heights = c(1.2, 4))
dev.off()
