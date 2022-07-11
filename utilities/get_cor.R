library(data.table, quietly=TRUE)
library(dplyr, quietly=TRUE)
library(boot)
library(stats)
library(parallel)
library(tidyverse)
library(doMC)

doMC::registerDoMC(cores = 17)
exp_id_path = "/home/dhthutrang/ENCODE/utilities/exp_id.txt"
dexseqcount_path = "/home/dhthutrang/ENCODE/mRNA_seq/dexseqcount/res"
all_pairs.exp_path = "/home/dhthutrang/ENCODE/flank/110722_1/all_pairs.exp.RDS"
filter_genes.exp_path = "/home/dhthutrang/ENCODE/utilities/combined_df_exon_90_final.RDS"
all_pairs.exp_flt_path = "/home/dhthutrang/ENCODE/flank/110722_1/all_pairs.exp_flt_90.RDS"

all_pairs.his_path = "/home/dhthutrang/ENCODE/flank/110722_1/all_pairs.his_list.RDS"
histone_type_list = c("H3K27ac", "H3K27me3", "H3K36me3", "H3K4me3", "H3K9me3")
filter_genes.his_path = "/home/dhthutrang/ENCODE/utilities/combined_df_exon_90_final.RDS"
all_pairs.his_flt_path = "/home/dhthutrang/ENCODE/flank/110722_1/all_pairs.his_list_flt_90_manorm.RDS"

all_res_list.pearcor_p_path = "/home/dhthutrang/ENCODE/flank/110722_1/all_res_list.pearcor_p_90_manorm.RDS"
all_res_list.pearcor_r_path = "/home/dhthutrang/ENCODE/flank/110722_1/all_res_list.pearcor_r_90_manorm.RDS"

get_colname <- function(filename_list, option='his'){
  name = sapply(filename_list, function(x) strsplit(x, split='/'))
  name = sapply(name, function(x) x[length(x)][[1]])
  if (option=="his"){
    name = sapply(name, function(x) strsplit(x, split="\\.")[[1]][1])
  }
  else if (option=="exp"){
    name = sapply(name, function(x) strsplit(x, split=".res.csv")[[1]][1])
  }
  name = sapply(name, function(x) paste(sort(strsplit(x, split="_")[[1]] ), collapse = "_"))
  names(name) = NULL
  return(name)
}

#= PREPARE EXP FILE (1 FOR ALL HIS TYPES) ====
print("=============== PREPARE EXP FILE (1 FOR ALL HIS TYPES) ===============")
all_pairs.exp = list.files(dexseqcount_path, full.names = TRUE, pattern='_res.csv')
get_all_pairs.exp <- function(all_pairs.exp){
  colname_exp = c("gene_id", "exon_id", get_colname(all_pairs.exp, "exp"))
  pair.exp_list = vector("list", length(all_pairs.exp))
  for (i in 1:length(all_pairs.exp)){
    print(paste("Pair: ", i, ': ', all_pairs.exp[[i]], sep=''))
    pair.exp = fread(all_pairs.exp[[i]])
    #fwrite(exp_id, "/home/dhthutrang/ENCODE/utilities/exp_id.txt", col.names = FALSE, row.names = FALSE, quote = FALSE, sep='\t')
    pair.exp_list[[i]] = pair.exp[,c("stat", "padj")]
  }
  pair.exp_list = lapply(pair.exp_list,
                         function(x) {
                           x = x %>%
                             mutate(
                               exp = dplyr::if_else(padj <= 0.05 & !is.na(padj),
                                                    true = stat, false = 0.0)) %>%
                             dplyr::select(exp)
                         })
  pair.exp_list = as.data.frame(pair.exp_list)
  exp_id = fread(exp_id_path, sep = '\t', quote=FALSE, header = FALSE)
  pair.exp_list = as.data.frame(cbind(exp_id, pair.exp_list))
  pair.exp_list = pair.exp_list[order(pair.exp_list$V1), ]
  colnames(pair.exp_list) = colname_exp
  return(as.data.table(pair.exp_list))
}


# FILTER WITH NEW CORRECTED GENES AND FIRST EXON
filter_genes = function(df, filter_genes_path="combined_df_exon.RDS", filter="deu"){
  combined_df_exon = readRDS(filter_genes_path)
  corrected_genes = unique(unlist(combined_df_exon[combined_df_exon['deu'] == T, "gene_id"]))

  filtered_df = as.data.frame(df[df$gene_id %in% corrected_genes, ]) #Filter by corrected genes[only 10 percs]
  combined_df = as.data.frame(combined_df_exon[combined_df_exon$gene_id %in% corrected_genes, ])
  filtered_df = filtered_df[order(filtered_df$gene_id), ] #Order so rownames overlap
  combined_df = combined_df[order(combined_df$gene_id), ]
  filtered_df = filtered_df[combined_df$first_exon == F, ] #Filter out first exons
  combined_df = combined_df[combined_df$first_exon == F, ]

  final_filtered_df = filtered_df #Filter by deu/dhm
  if (filter == "deu") final_filtered_df[unlist(combined_df[[filter]] == F), 3:length(final_filtered_df)] = 0

  #------Check if filter by DEU work------
  print("========= CHECK FILTER =========")
  print(paste("Gene IDs match:", table(final_filtered_df$gene_id == filtered_df$gene_id)))
  print(paste("Overlap before after:", table(final_filtered_df == filtered_df)))
  return(as.data.table(final_filtered_df))
}

# # === EXECUTE get_all_pairs.exp ======
# all_pairs.exp = get_all_pairs.exp(all_pairs.exp)
# saveRDS(all_pairs.exp, all_pairs.exp_path)
# all_pairs.exp_ = readRDS(all_pairs.exp_path)

# # === EXECUTE FILTER get_all_pairs.exp ======
# all_pairs.exp_flt_90 = filter_genes(all_pairs.exp_, filter_genes_path=filter_genes.exp_path, filter="deu")
# saveRDS(all_pairs.exp_flt_90, all_pairs.exp_flt_path)


all_pairs.exp_flt_90 = readRDS(all_pairs.exp_flt_path)

# Frequency of DEU exons in each pairwise comp
# x = apply(all_pairs.exp[, 3:ncol(all_pairs.exp)], 2, function(x) as.data.frame(table(table(all_pairs.exp$gene_id[x > 0]))))
# y = do.call('rbind', x) %>%
#   dplyr::group_by(Var1) %>%
#   dplyr::summarise(n = sum(Freq)) %>%
#   ungroup() %>%
#   mutate(perc = round(n/sum(n), 5)*100)
# head(y)
# ggplot(y, aes(x=Var1, y=n)) +
#   geom_col()

# = PREPARE HIS FILE (6 TOTAL) ====
print("=============== PREPARE HIS FILE (6 TOTAL) ===============")
get_all_pairs.his <- function(all_pairs.his, his){
  pair.his_list = vector("list", length(all_pairs.his))
  for (i in 1:length(all_pairs.his)){
    # for (i in 1:1){
    print(paste("Pair: ", i, sep=''))
    pair.his = fread(all_pairs.his[[i]])
    id = as.data.frame(do.call(rbind, lapply(pair.his$V9, function(x) strsplit(x, split='"', fixed=T)[[1]][c(2, 6)])))
    colnames(id) = c('gene', 'exon')
    pair.his = pair.his %>%
      dplyr::mutate(
        gene = id$gene, exon = id$exon, type = V3,
        p_val = as.numeric(as.character(V11)),
        m_val = dplyr::if_else(
          p_val <= 0.05, 
          true = as.numeric(as.character(V10)), false = 0
          ) # Transform values to absolute
      ) 
      pair.his %>%
      dplyr::select(gene, exon, m_val) %>%
      dplyr::group_by(gene, exon) %>% 
      dplyr::summarise_all(max, na.rm=T) %>%
      dplyr::na_if(., -Inf) %>% 
      dplyr::na_if(., Inf) %>%
      dplyr::ungroup() 
    # print(head(pair.his))
    
    if (i == 1) pair.his_id = pair.his[, c('gene', 'exon')]
    pair.his = pair.his %>% dplyr::select(-gene, -exon)
    pair.his_list[[i]] = pair.his
  }
  lapply(pair.his_list, function(x) print(dim(x)))
  pair.his_list = as.data.frame(cbind(pair.his_id, as.data.frame(do.call(cbind, pair.his_list))))
  # print(dim(pair.his_list))
  # print(head(pair.his_list))
  saveRDS(pair.his_list, paste('pair.his_list_', his, '.RDS', sep=''))
  return(pair.his_list)
}

get_all_pairs.his_list <- function(histone_type_list){
  all_pairs.his_list = vector("list", length(histone_type_list))
  for (j in 1:length(histone_type_list)){
    his = histone_type_list[[j]]
    all_pairs.his = list.files(paste("/home/dhthutrang/ENCODE/chip_seq", his, "flank/fl", sep='/'), pattern = '.txt', full.names = TRUE)
    # print(all_pairs.his)
    colname_his = c("gene_id", "exon_id", get_colname(all_pairs.his, "his")) 
    all_pairs.his.sig = get_all_pairs.his(all_pairs.his, his)
    colnames(all_pairs.his.sig) = colname_his
    
    all_pairs.his_list[[j]] = as.data.table(all_pairs.his.sig)
  }
  return(all_pairs.his_list)
}

filter_all_his_list <- function(his_list, histone_type_list, filter_genes_path){
  all_filtered_df = vector("list")
  for (i in 1:length(histone_type_list)){
    histone = histone_type_list[i]
    his_df = his_list[[i]]
    print(i)
    print(histone)
    all_filtered_df[[i]] = filter_genes(df = his_df, filter_genes_path = filter_genes_path, filter = histone)
  }
  names(all_filtered_df) = histone_type_list
  return(all_filtered_df)
}

# === EXECUTE get_all_pairs.his_list ======
all_pairs.his_list = get_all_pairs.his_list(histone_type_list)
saveRDS(all_pairs.his_list, all_pairs.his_path)
all_pairs.his_list_ = readRDS(all_pairs.his_path)
names(all_pairs.his_list_) = histone_type_list

# === EXECUTE file get_all_pairs.his_list ======
all_pairs.his_list_flt_90 = filter_all_his_list(all_pairs.his_list_, histone_type_list, filter_genes.his_path)
saveRDS(all_pairs.his_list_flt_90, all_pairs.his_flt_path)


all_pairs.his_list_flt_90 = readRDS(all_pairs.his_flt_path)

#== Plotting genes ----
get_plot_dfs = function(gene, tissue_pair){
  return(
    lapply(all_pairs.his_list_flt_90, function(x)
      df = data.frame(
        DEU=all_pairs.exp_flt_90[[tissue_pair]][all_pairs.exp_flt_90$gene_id == gene], 
        DHM=x[[tissue_pair]][x$gene_id == gene]) 
    )
  )
}
# LMNB1_dfs = get_plot_dfs("LMNB1", "neuronalstemcell_pancreas")
# FGFR2_dfs = get_plot_dfs("FGFR2", "mesenchymalstemcell_sigmoidcolon")

plot_cor = function(dfs, title_prefix){
  name_dfs = names(dfs)
  plots = list()
  for (i in 1:length(dfs)){
    df = dfs[[i]]
    plots[[i]] = ggscatter(data=df, x = 'DEU', y = 'DHM',
                           add = "reg.line", 
                           add.params = list(color = "#5D2ED2", fill = "lightgray"),
                           conf.int = TRUE, 
                           font.label = c(8, 'plain', 'red'),
                           title = paste(title_prefix, name_dfs[i], sep='\n')) +
      stat_cor(method = "pearson", label.x=min(df$DEU) + 0.5, label.y = max(df$DHM) + 0.5) +
      theme(aspect.ratio = 0.8,
            title = element_text(size=9),
            axis.title = element_text(face = 'bold', size=12))
  }
  names(plots) = name_dfs
  return(plots)
}
# LMNB1_plots = plot_cor(LMNB1_dfs, "LMNB1\nNeuronal stem cell - Pancreas")
# FGFR2_plots = plot_cor(FGFR2_dfs, "FGFR2\nMesenchymal stem cell - Sigmoid colon")
# LMNB1_plots$H3K27ac
# plot_list = list(LMNB1_plots[[1]], LMNB1_plots[[3]], LMNB1_plots[[4]], FGFR2_plots[[2]])
# png("temp.png", width=7, height = 7, units = 'in', res=300)
# ggarrange(plotlist = plot_list, nrow = 2, ncol=2, align = 'hv')

#== Get all_pairs.his_list_flt binary ----
his_to_binary <- function(all_pairs.his_list_flt){
  print(length(all_pairs.his_list_flt))
  all_pairs.his_list_flt_bin = list()
  for (i in 1:length(all_pairs.his_list_flt)){
    print(paste("CHECK ", i))
    dhm = all_pairs.his_list_flt[[i]][, 3:ncol(all_pairs.his_list_flt[[i]])]
    dhm = dhm > 0
    print(head(dhm))
    dhm_bin = apply(dhm, 1, function(x) Reduce(function(a,b) a|b, na.omit(x)))
    all_pairs.his_list_flt_bin[[i]] = dhm_bin
  }
  lapply(all_pairs.his_list_flt_bin, function(x) print(dim(x)))
  all_pairs.his_list_flt_bin = cbind(all_pairs.his_list_flt[[1]][,1:2], do.call(cbind, all_pairs.his_list_flt_bin))
  return(all_pairs.his_list_flt_bin)
}
# all_pairs.his_list_flt_90_bin = his_to_binary(all_pairs.his_list_flt_90)
# saveRDS(all_pairs.his_list_flt_90_bin, "/home/dhthutrang/ENCODE/utilities/all_pairs.his_list_flt_90_bin.RDS")

# #===== CORRELATION WITH RANDOMIZATION =====
# ------------ Execute analysis ------------
all_genes = fread("gene_id.txt", header = FALSE) #TODO 
all_genes = unique(all_pairs.exp_flt_90$gene_id)
p_value_calculator <- function(r, nrow){
  P <- r*sqrt(nrow-2)/sqrt(1-r*r)
  P <- 2*pt(-abs(P), nrow-2)
  return(P)
}
pearcor_p <- function(exp, his, cor_type = "pearson"){
  if (length(unique(exp)) > 1 & length(unique(his)) > 1){
    p_val = p_value_calculator(cor(exp, his, method = cor_type), nrow = length(exp))
    # p_val = p.adjust(p_val, method = "fdr", n=ncol(all_genes))
    return(p_val)
  }
  else {
    return(NA)
  }
}

pearcor_r <- function(exp, his, n_points, cor_type = "pearson"){
  df = as.data.frame(cbind(exp, his))
  n_sep_point = nrow(unique(df))
  if (0 %in% apply(df, 1, unique)) n_sep_point = n_sep_point - 1
  if (n_sep_point >= 3 & length(unique(exp)) > 1 & length(unique(his)) >= n_points){
    r_val = cor(exp, his, method = cor_type)
    return(r_val)
  }
  else {
    return(NA)
  }
}

analyze_array <- function(all_pairs.exp, all_pairs.his, option = "p", n_points, all_genes, cor = "pearson"){
  all_res_pair = vector("list", ncol(all_pairs.exp) - 2)
  subset_name = colnames(all_pairs.his)
  # print(head(subset_name))
  # print(head(all_pairs.exp))
  colnames(all_pairs.exp) = gsub('trophoblastcell', 'trophoblast', colnames(all_pairs.exp))
  all_pairs.exp_subset = all_pairs.exp[, ..subset_name]

  #for (i in 1:n_pairs){
  all_res_pair <- foreach( i=1:(length(subset_name)-2), .combine='c', .packages=c('dplyr') ) %dopar% {
    print(paste("Pair: ", i, sep=''))
    exp = all_pairs.exp_subset[[i+2]]
    his = all_pairs.his[[i+2]]
    data_table = as.data.table(cbind(exp, his))
    # head(data_table)

    if (option == "p"){
      res_table = data_table %>%
        dplyr::group_by(all_pairs.exp$gene_id) %>%
        dplyr::summarise(res = pearcor_p(exp, his, cor_type = cor)) %>%
        dplyr::select(res)
    }
    else if (option == "r") {
      res_table = data_table %>%
        dplyr::group_by(all_pairs.exp$gene_id) %>%
        dplyr::summarise(res = pearcor_r(exp, his, n_points, cor_type = cor)) %>%
        dplyr::select(res)
    }
    #all_res_pair[[i]] = res_table
  }
  all_res_pair = as.data.table(all_res_pair)
  all_res_pair = cbind(all_genes, all_res_pair)
  colnames(all_res_pair) = c('gene_id', subset_name[3:length(subset_name)])
  # print(head(all_res_pair))
  return(as.data.frame(all_res_pair))
}
analyze_array_list <- function(all_pairs.exp, all_pairs.his_list, method = "p", n_points=2, cor = "pearson"){
  all_res_list = vector("list", length(histone_type_list)-1 )
  all_genes = unique(all_pairs.exp$gene_id)
  for (j in 1:length(histone_type_list)){
    print(paste("Histone: ", histone_type_list[[j]], sep = ''))
    all_pairs.his = all_pairs.his_list[[j]]
    if (method == "p") {
      all_res_pair = analyze_array(
        all_pairs.exp = all_pairs.exp, 
        all_pairs.his = all_pairs.his, 
        option = "p", 
        all_genes=all_genes, 
        cor = cor)
    }
    else if (method == "r") {
      print("CHECK 0")
      all_res_pair = analyze_array(
        all_pairs.exp = all_pairs.exp, 
        all_pairs.his = all_pairs.his, 
        option = "r", 
        n_points=n_points, 
        all_genes=all_genes, 
        cor = cor)
    }
    all_res_list[[j]] = all_res_pair
  }
  return(all_res_list)
}

print("Pearsons-p correlation")
all_res_list.pearcor_p = analyze_array_list(
  all_pairs.exp = all_pairs.exp_flt_90,
  all_pairs.his_list = all_pairs.his_list_flt_90, 
  method = "p",
  cor = "pearson")
saveRDS(all_res_list.pearcor_p, all_res_list.pearcor_p_path)
# 
all_res_list.pearcor_r = analyze_array_list(
  all_pairs.exp = all_pairs.exp_flt_90, 
  all_pairs.his_list = all_pairs.his_list_flt_90, 
  method = "r", 
  cor = "pearson"
  )
saveRDS(all_res_list.pearcor_r, all_res_list.pearcor_r_path)
