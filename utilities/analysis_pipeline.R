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
##
## ---------------------------

library("optparse", quietly = TRUE)

# ==== PREPARATION ====
# Parse arguments for input/output settings
option_list <- list(
  make_option("--general_analysis_results",
    type = "character",
    help = "path to the folder where the general results from the analysis are stored and shared between processes",
    metavar = "character",
    default = "general_analysis_results"
  ),
  make_option(c("-f", "--datasets_path"),
    type = "character",
    help = "path to the folder where the general results from the analysis are stored for all used dataset",
    metavar = "character", 
    default = "general_analysis_results/datasets"
  ),
  make_option("--multi_datasets_results",
    type = "character",
    help = "path to the folder where the results from the comparing different datasets are stored",
    metavar = "character",
    default = "multi_datasets_results"
  ),
  make_option("--epigenome_information",
    type = "character",
    help = "path to a ; separated dataframe with Epigenome, Potency, Type, Origin, Life stage, Official name are",
    metavar = "character",
    default = "general_analysis_results/epigenome_info.csv"
  ),
  make_option("--analysis_methods",
    type = "character",
    help = "path to the script containing analysing methods",
    metavar = "character",
    default = "utilities/analysis_methods.R"
  ),
  make_option("--analysis_methods_GO",
    type = "character",
    help = "path to the script containing analysing methods for Gene Ontology Enrichment Analysis",
    metavar = "character",
    default = "utilities/analysis_methods_GO.R"
  )
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

general_analysis_results <- opt$general_analysis_results
dataset_paths <- dir(opt$datasets_path, full.names = TRUE)
sapply(dataset_paths, function(path_) dir.create(file.path(path_, "res"), showWarnings = FALSE))

multi_datasets_results_path <- opt$multi_datasets_results_path
dir.create(file.path(multi_datasets_results_path), showWarnings = FALSE)

source(opt$analysis_methods)
source(opt$analysis_methods_GO)
shinyGO_result_path <- paste(general_analysis_results, "shinyGO_epispliced_genes.csv", sep = "/")
universe_genes_path <- paste(general_analysis_results, "all_DEU_genes.RDS", sep = "/")

histone_type_list <- c("H3K27ac", "H3K27me3", "H3K36me3", "H3K4me3", "H3K9me3")
epigenomes_annot <- read.csv(opt$epigenome_information, header = TRUE, sep = ";", row.names = 1)
official_name <- epigenomes_annot$Official.name


#==== ANALYSIS PIPELINE ====
run_pipeline <- function(dataset_path, dataset_index, r_sig = 0.5) {
  #' This function executes the analysis pipeline with the helping functions defined in "METHODS TO COLLECT TISSUE-
  #' SPECIFIC EPISPLICED GENES"
  #' 
  #' @param dataset_path a list of paths to the results of get_cor.R. Each path leads to a folder where 4 dataframes
  #' containing DEU statistics, DHM M-values, correlation p-values and correlation r-values are stored in RDS format
  #' @param dataset_index the index of the dataset_path being analyse. Used for reporting during parallel process only
  #' @param r_sig the threshold for correlation R-values to define epispliced genes
  #' @return a list containg the analysed results for each dataset. Each result contains 2 lists: 'info' where paths 
  #' and parameters set for the analysis are stored, and 'results' where the generated lists and dataframes are stored
  #'  
  
  info <- list()
  results <- list()
  dataset_name <- sapply(dataset_paths, function(path_) {
    name <- strsplit(path_, "/")[[1]]
    return(name[length(name)])
    })
  absolute_values <- strsplit(dataset_name, "_")[[1]][1]
  correlation_type <- strsplit(dataset_name, "_")[[1]][2]

  info["absolute_values"] <- absolute_values
  info["correlation_type"] <- correlation_type

  print(paste(
    "Currently processing: ", dataset_path, " with options: Absolute values - ", absolute_values,
    " Correlation type - ", correlation_type,
    sep = ""
  ))


  # Define paths
  print("..... Define paths")
  dataset_path <- dataset_path
  res_p_value_path <- paste(dataset_path, "p_val_dfs.RDS", sep = "/")
  res_r_value_path <- paste(dataset_path, "r_val_dfs.RDS", sep = "/")
  DEU_path <- paste(dataset_path, "DEU_df.RDS", sep = "/")
  DHM_path <- paste(dataset_path, "DHM_dfs.RDS", sep = "/")

  info["dataset_path"] <- dataset_path
  info["res_p_value_path"] <- res_p_value_path
  info["res_r_value_path"] <- res_r_value_path
  info["DEU_path"] <- DEU_path
  info["DHM_path"] <- DHM_path


  # Correct R
  print("..... Correct R")
  all_p <- readRDS(res_p_value_path)
  all_padj <- get_p_adj(p_val_list = all_p)

  results["all_padj"] <- list(all_padj)


  # Filter R values by p-values
  print("..... Filter R values by p-values")
  all_r <- readRDS(res_r_value_path)
  all_r_adjusted <- filter_by_p(r_df_list = all_r, p_df_list = all_padj)

  results["all_r_adjusted"] <- list(all_r_adjusted)


  # Get data as controls
  print("..... Get data as controls")
  DEU_df <- readRDS(DEU_path)
  all_r_adjusted_DEU <- DEU_df %>%
    dplyr::select(-exon_id) %>%
    dplyr::group_by(gene_id) %>%
    dplyr::mutate_all(abs) %>%
    dplyr::summarise_all(sum) %>%
    as.data.frame()
  all_r_adjusted_DEU <- list(all_r_adjusted_DEU)

  DHM_df <- readRDS(DHM_path)
  all_r_adjusted_DHM <- lapply(DHM_df, function(x) {
    return(x %>% dplyr::select(-exon_id) %>%
      dplyr::group_by(gene_id) %>%
      dplyr::mutate_all(abs) %>%
      dplyr::summarise_all(sum) %>%
      as.data.frame())
  })
  names(all_r_adjusted_DHM) <- histone_type_list

  results["all_r_adjusted_DEU"] <- list(all_r_adjusted_DEU)
  results["all_r_adjusted_DHM"] <- list(all_r_adjusted_DHM)


  # Summarize significant genes for each tissue pair
  print("..... Summarize significant genes for each tissue pair")
  all_padj_sig <- get_sig_gene_list(r_df_list = all_r_adjusted, method = "pearcor_r", r_sig = r_sig)
  all_padj_sig_DEU <- get_sig_gene_list(r_df_list = all_r_adjusted_DEU, method = "DEU", DEU_sig = 2)
  all_padj_sig_DHM <- get_sig_gene_list(r_df_list = all_r_adjusted_DHM, method = "DHM", DHM_sig = 1)

  results["all_padj_sig"] <- list(all_padj_sig)
  results["all_padj_sig_DEU"] <- list(all_padj_sig_DEU)
  results["all_padj_sig_DHM"] <- list(all_padj_sig_DHM)
  

  # Get significant genes for each tissue 
  print("..... Get significant genes for each tissue")
  tissue_type_list <- get_tissue_specific_indices(r_sig_list = all_padj_sig)
  tissue_type_list_DEU <- tissue_type_list[1:1]
  names(tissue_type_list_DEU) <- "DEU" 

  TSEGs <- get_TSEGs(r_sig_list = all_padj_sig, tissue_spec_ind = tissue_type_list)
  DEU_genes <- get_TSEGs(r_sig_list = all_padj_sig_DEU, tissue_spec_ind = tissue_type_list_DEU)
  DHM_genes <- get_TSEGs(r_sig_list = all_padj_sig_DHM, tissue_spec_ind = tissue_type_list)

  all_TSEGs <- unique(unlist(unlist(TSEGs)))
  print(paste("..... ..... Number of total significant genes: ", length(all_TSEGs), sep = ""))

  results["TSEGs"] <- list(TSEGs)
  results["DEU_genes"] <- list(DEU_genes)
  results["DHM_genes"] <- list(DHM_genes)
  results["all_TSEGs"] <- list(all_TSEGs)
  

  # If there is pos and neg
  print("..... Extra: Gather postive and negative R-values")
  if (absolute_values == "True"){
    all_padj_sig_pos = get_sig_gene_list(r_df_list = all_r_adjusted, method = "pearcor_r", 
                                         r_sig = 0, compare_threshold = "greater")
    all_padj_sig_neg = get_sig_gene_list(r_df_list = all_r_adjusted, method = "pearcor_r", 
                                         r_sig = 0, compare_threshold = "lesser")
  
    TSEGs_pos = get_TSEGs(r_sig_list = all_padj_sig_pos, tissue_spec_ind = tissue_type_list)
    TSEGs_neg = get_TSEGs(r_sig_list = all_padj_sig_neg, tissue_spec_ind = tissue_type_list)

    results['all_padj_sig_pos'] = list(all_padj_sig_pos)
    results['all_padj_sig_neg'] = list(all_padj_sig_neg)
    results['TSEGs_pos'] = list(TSEGs_pos)
    results['TSEGs_neg'] = list(TSEGs_neg)
  }
  else {
    results['all_padj_sig_pos'] = NULL
    results['all_padj_sig_neg'] = NULL
    results['TSEGs_pos'] = NULL
    results['TSEGs_neg'] = NULL
  }

  # Return outcome with results stored in dataframes and lists
  outputs <- list(info, results)
  names(outputs) <- c("info", "results")
  return(outputs)  
}

post_analysis <- function(analysis_results_){
  #' This function generates dataframes and lists using functions defined in "METHODS TO PREPARE PLOTS & TABLES DATA"
  #' to prepare for figures/tables generation 
  #' 
  #' @param analysis_results_ a list containg the analysed results for each dataset. Each result contains 2 lists: 
  #' 'info' - where paths and parameters set for the analysis are stored, and 'results' - where the generated lists and
  #' dataframes are stored.
  #' @return a list containg further analysed results for each dataset
  #' 
  
  info <- analysis_results_[["info"]]
  results <- analysis_results_[["results"]]
  dataset_path <- info[["dataset_path"]]

  #---- Get cumulative plot for R-values (Figure 3A) ----
  print("..... Get data for cumulative plot of R-values (Figure 3A)")
  # Prepare datasets
  all_r <- readRDS(info[["res_r_value_path"]])
  all_r_cor <- prep_plot_ecdf(r_df_list = all_r)
  all_r_cor$type <- "padj > 0.05"
  
  all_r_adjusted <- results[["all_r_adjusted"]] 
  all_r_adjusted_cor <- prep_plot_ecdf(r_df_list = all_r_adjusted)
  all_r_adjusted_cor$type <- "padj <= 0.05"
  all_r_compare <- rbind(all_r_cor, all_r_adjusted_cor)

  results["all_r_compare"] <- list(all_r_compare)

  # Plot figure
  plot_ecdf(all_r_compare = all_r_compare, dataset_path = dataset_path)


  #---- Get and report the number of epispliced genes for each tissue (Table 2) ----
  print("..... Report the number of epispliced genes for each tissue (Table 2)")
  # Prepare datasets
  TSEGs <- results[["TSEGs"]]
  no_sig_gene_by_tissue <- get_TSEGs_no(TSEGs_list = TSEGs)
  no_sig_gene_by_tissue$`Total epispliced genes (with overlaps)` <- 
    rowSums(apply(no_sig_gene_by_tissue[, 2:6], 2, as.numeric), na.rm = TRUE)
  no_sig_gene_by_tissue$`Total epispliced genes (without overlaps)` <- 
    sapply(get_nonoverlap_TSEG_set(TSEGs), length)
  no_sig_gene_by_tissue <- as.data.frame(no_sig_gene_by_tissue)

  results["no_sig_gene_by_tissue"] <- list(no_sig_gene_by_tissue)

  # Plot figure
  export_no_sig_gene_by_tissue(no_sig_gene_by_tissue = no_sig_gene_by_tissue, dataset_path = dataset_path)


  #---- Get data for the heatmap representing the distance between tissue pairs (Figure 5 and Subfigure S3) ----
  print("..... Get tissue pairwise distance (Figure 5 & Subfigure S3)")
  # Prepare datasets
  DEU_genes <- results[["DEU_genes"]]
  DHM_genes <- results[["DHM_genes"]] 

  sim_mat <- get_similarity_matrix(TSEGs, histone_type_list)
  sim_mat_DEU <- get_similarity_matrix(DEU_genes, "DEU")
  sim_mat_DHM <- get_similarity_matrix(DHM_genes, histone_type_list)
  sim_mat_DEU_DHM <- c(sim_mat_DEU, sim_mat_DHM)
  
  results["sim_mat"] <- list(sim_mat)
  results["sim_mat_DEU_DHM"] <- list(sim_mat_DEU_DHM)

  # Plot figure
  plot_heatmap(sim_mat = sim_mat, sim_mat_DEU_DHM = sim_mat_DEU_DHM, dataset_path = dataset_path)


  #---- Plot the DEU statistics, DHM M-values and correlation for example genes (Figure 3B) ----
  DEU_df <- readRDS(info[["DEU_path"]])
  DHM_dfs <- readRDS(info[["DHM_path"]])
  LMNB1_dfs <- get_plot_dfs(DEU_df = DEU_df, DHM_dfs = DHM_dfs, gene = "LMNB1", 
                            tissue_pair = "neuronalstemcell_pancreas")
  FGFR2_dfs <- get_plot_dfs(DEU_df = DEU_df, DHM_dfs = DHM_dfs, gene = "FGFR2",
                            tissue_pair = "mesenchymalstemcell_sigmoidcolon")

  cor_plots <- list()
  LMNB1_cor_plots <- plot_cor(dfs = LMNB1_dfs,
                              plot_title = paste("LMNB1", "Neuronal stem cell - Pancreas"))
  FGFR2_cor_plots <- plot_cor(dfs = FGFR2_dfs,
                              plot_title = paste("FGFR2", "Mesenchymal stem cell - Sigmoid colon"))
    
  # Plot DEU, DHM and correlation for LMNB1 (only H3K27ac, H3K36me3, H3K4me3) and FGFR2 (only H3K27me3)
  plot_list <- list(LMNB1_cor_plots[[1]], LMNB1_cor_plots[[3]], LMNB1_cor_plots[[4]], FGFR2_cor_plots[[2]])
  tiff(paste(dataset_path, "res", "example_genes.tiff", sep = "/"), width = 7, height = 4, units = "in", res = 300)
  print(ggarrange(plotlist = plot_list, ncol = 2, nrow = 2, align = "hv"))
  dev.off()

  # ---- Return outcome with updated resulting dataframes and lists ----
  outputs <- list(info, results)
  names(outputs) <- c("info", "results")
  return(outputs)
}


GO_term_enrichment_analysis <- function(analysis_results_, universe_genes_path, shinyGO_result_path){
  #' This function generates GO Terms Enrichment Analysis (GOEA) results using functions defined in "METHODS TO 
  #' PERFORM GENE ONTOLOGY ENRICHMENT ANALYSIS (GOEA)" in GO_analysis_methods.R
  #' 
  #' @param analysis_results_ a list containg the analysed results for each dataset. Each result contains 2 lists: 
  #' 'info' - where paths and parameters set for the analysis are stored, and 'results' - where the generated lists and
  #' dataframes are stored.
  #' @return a list containg further analysed results for each dataset
  #' 
  
  #---- Prepare data for the gene annotations enrichment analysis ----
  print("..... Get data for GO terms enrichment analysis")
  info <- analysis_results_[["info"]]
  results <- analysis_results_[["results"]]
  dataset_path <- info[["dataset_path"]]

  universe_genes <- readRDS(universe_genes_path) 
  universe_genes <- unique(unlist(lapply(universe_genes, function(x) strsplit(x, ";")[[1]][1])))
  universe_genes_entrez <- get_entrez_id(universe_genes)
  universe_genes_entrez <- universe_genes_entrez[!is.na(universe_genes_entrez)]

 # !!! IMPORTANT: Might run into "database disk image is malformed" error. Use biomartCacheClear() to clear cache !!!
  signs <- c("none", "neg", "pos")
  GO_result_list <- list()
  for (i in seq(length(signs))){
    GO_result <- tryCatch({
      run_GO_analysis(results, universe_genes_entrez = universe_genes_entrez, sign = signs[i])
      }, 
      error = function(e) {print("Check analyse_multiple_GO:run_GO_analysis()")})
    GO_result_list[[i]] <- GO_result
  }
  names(GO_result_list) <- signs 
  results[["GO_result_list"]] <- GO_result_list


  #---- Plot the results for gene annotations enrichment analysis (Figure 6) ----
  print("..... Plot enrichment from negative and postive correlation (Subfigre S4)")
  get_dotplot_signed_results(GO_result_list = GO_result_list, 
    filename = paste(dataset_path, "res", "GO_enrichment_signed.tiff", sep = "/"))

  print("..... Plot enrichment from main results (Figure 6)")
  get_dotplot(GO_result_list = GO_result_list, gene_ratio = 0.01, 
    filename = paste(dataset_path, "res", "GO_enrichment.tiff", sep = "/"))

  print("..... Plot enrichment from main results for all epispliced genes (Figure 6)")
  shinyGO_result <- fread(shinyGO_result_path, header = TRUE)
  get_dotplot_all_epispliced_genes(annot_df = shinyGO_result, 
    filename = paste(dataset_path, "res", "GO_enrichment_all_epispliced_genes.tiff", sep = "/"))
  
  # ---- Return outcome with updated resulting dataframes and lists ----
  outputs <- list(info, results)
  names(outputs) <- c("info", "results")
  return(outputs)
}

#==== EXECUTE PIPELINE ====
res_all_datasets <- foreach(i = seq(length(dataset_paths)), 
  .combine = "list", .multicombine = TRUE, .errorhandling = 'pass') %dopar% { 
  output <- run_pipeline(
    dataset_path = dataset_paths[i], 
    dataset_index = i
    )
}
saveRDS(res_all_datasets, paste(multi_datasets_results_path, "res_all_datasets.RDS", sep = "/"))

res_all_datasets_extended <- foreach(i = seq(length(res_all_datasets)), 
  .combine = "list", .multicombine = TRUE, .errorhandling = 'pass') %dopar% { 
  output <- post_analysis(analysis_results_ = res_all_datasets[[i]])
}
saveRDS(res_all_datasets_extended, paste(multi_datasets_results_path, "res_all_datasets_extended.RDS", sep = "/"))

# Gene enrichment analysis must be called sequentially because of cache storage size limit
for (result in res_all_datasets_extended){
  GO_term_enrichment_analysis(
    analysis_results_ = result, 
    universe_genes_path = universe_genes_path, 
    shinyGO_result_path = shinyGO_result_path
    )
}

#==== METHODS TO COMPARE THE RESULTS FROM DIFFERENT DATASETS ====
library("VennDiagram")
draw_venn <- function(gene_sets){
  #' Draw the Venn diagram for the epipliced genes overlap between different datasets in Subfigure S1B
  #' 
  #' @param gene_sets a named list of epispliced genes from each dataset 
  #' @return None
  #' 
  
  n_groups <- length(gene_sets)
  venn_diag = venn.diagram(
          x = gene_sets, category.names = names(gene_sets),
          filename = paste(multi_datasets_results_path, "compare_datasets_venn.tiff", sep = "/"), output=TRUE,
          imagetype="tiff" , height = 5 , width = 5 , units = "in", resolution = 300,
          fontfamily = "sans", print.mode = "percent", sigdigs = 2,
          fill = brewer.pal(length(gene_sets), "Pastel2")[seq(n_groups)],
          cat.cex = 0.8, cat.fontface = "bold", cat.default.pos = "outer", cat.fontfamily = "sans", 
          cat.pos = rep(0, n_groups), lwd = 2, lty = 'blank', cex = 0.8
  )
}

get_name_list <- function(analysis_results_){
  #' Extract names of datasets being compared against each other into a list 
  #' 
  #' Given a list of results for different datasets, get their names which should follow the format 
  #' "absolute_values"_"correlation_type". Absolute_values is either Abs (use absolute values of DEU/DHM values) or 
  #' True (use true values) and correlation_type is either Pearson or Spearman. 
  #' @param analysis_results_ a list containg the analysed results for each dataset. Each result contains 2 lists: 
  #' 'info' - where paths and parameters set for the analysis are stored, and 'results' - where the generated lists and
  #' dataframes are stored.
  #' 
  
  name_list_ <- lapply(analysis_results_, function(x) 
    paste(x[["info"]][["absolute_values"]], x[["info"]][["correlation_type"]], sep = "_")
  )
  return(name_list_)
}

compare_all_TSEG <- function(analysis_results_){
  #' Compute the Jaccard distance between the list of epispliced genes identified for each dataset and produce the Venn
  #' diagram for the overlap in epispliced genes from different datasets in Subfigure S1B
  #' 
  #' @param analysis_results_ a list containg the analysed results for each dataset. Each result contains 2 lists: 
  #' 'info' - where paths and parameters set for the analysis are stored, and 'results' - where the generated lists and
  #' dataframes are stored.
  #' 
  
  # Gather all epispliced genes detected for each dataset 
  all_TSEG_list <- lapply(analysis_results_, function(x) x[["results"]][["all_TSEGs"]])
  name_list <- get_name_list(analysis_results_) 
  names(all_TSEG_list) = name_list

  # Compute the Jaccard distance between epispliced gene sets identified from different datasets
  distances <- as.data.frame(compute_jaccard(all_TSEG_list), row.names = name_list)
  colnames(distances) <- name_list
  print(distances)

  # Plot the Venn diagram to show the overlap in epispliced genes between datasets
  draw_venn(all_TSEG_list)
  
  return(distances)
}

compare_cors <- function(analysis_results_, cor_method){
  #' Examine the correlation between the R-values for DEU/DHM corrrelation computed for each dataset
  #' 
  #' @param analysis_results_ a list containg the analysed results for each dataset. Each result contains 2 lists: 
  #' 'info' - where paths and parameters set for the analysis are stored, and 'results' - where the generated lists and
  #' dataframes are stored.
  #' @param cor_method "pearson" or "spearman" correlation used to compare any two list of R-values
  #' 
  
  # Gather R-values from different datasets into a dataframe
  all_cors_list <- do.call("cbind", lapply(analysis_results_, function(x) melt(as.data.frame(x[['results']][['all_r_adjusted']][[1]]))$value))
  print(class(all_cors_list))
  print(head(all_cors_list))
  all_cors_list[is.na(all_cors_list)] = 0
  all_cors_list <- all_cors_list[apply(all_cors_list, 1, function(x) !all(x = 0)), ]
  name_list <- get_name_list(analysis_results_) 

  # Compute correlation
  print("Correlation between computed the set of DEU-DHM R-values from different datasets")
  cor_mat <- cor(all_cors_list, method = cor_method)
  colnames(cor_mat) <- name_list
  rownames(cor_mat) <- name_list 
  print(cor_mat)

  print("Correlation between computed the set of absolute values of DEU-DHM R-values from different datasets")
  cor_mat_abs <- cor(abs(all_cors_list), method = cor_method)
  colnames(cor_mat_abs) <- name_list
  rownames(cor_mat_abs) <- name_list
  print(cor_mat_abs)

  return(list(cor_mat, cor_mat_abs))
}

compare_all_TSEG(res_all_datasets_extended)
compare_cors(res_all_datasets_extended, "spearman")

