library(ggplot2)
library(reshape2)
library(ggpubr)
library(dplyr)
library(data.table)
setwd("/Users/dhthutrang/Documents/BIOINFO/Episplicing/ENCODE_episplicing/flank")
histone_type_list = c("H3K27ac", "H3K27me3", "H3K36me3", "H3K4me1", "H3K4me3", "H3K9me3")


all_pairs.exp = readRDS("all_pairs.exp.RDS")
all_pairs.exp = as.data.frame(all_pairs.exp)
all_pairs.his = readRDS("all_pairs.his_list.RDS")
all_pairs.his = lapply(all_pairs.his, as.data.frame)

pear_r = readRDS("all_res_list.pearcor_r.RDS")
pear_r_5 = readRDS("all_res_list.pearcor_r_5.RDS")
pear_r_10 = readRDS("all_res_list.pearcor_r_10.RDS")
pear_p_adj = readRDS("all_res_list.pearcor_padj.RDS")
pear_p = readRDS("all_res_list.pearcor_p.RDS")

filter_by_p <- function(pear_df, pear_p_df) {
  pear_df_filtered = list()
  for (i in 1:6){
    r_table = as.matrix(pear_df[[i]][, 2:ncol(pear_df[[i]])])
    p_table = as.matrix(pear_p_df[[i]][, 2:ncol(pear_p_df[[i]])])
    r_table[p_table > 0.05] = NA
    r_table = as.data.frame(r_table)
    r_table = cbind(pear_df[[i]][, 1:2], r_table)
    pear_df_filtered[[i]] = r_table
  }
  return(pear_df_filtered)
}
# pear_r_p = filter_by_p(pear_r, pear_p)
pear_r_p = filter_by_p(pear_r, pear_p_adj)
pear_r_5 = filter_by_p(pear_r_5, pear_p_adj)
pear_r_10 = filter_by_p(pear_r_10, pear_p_adj)

get_max_cor <- function(pear_df){
  all_index = sapply(pear_df, function(x) {
    idx = which(x == max(x[, 2:ncol(x)], na.rm = TRUE),  arr.ind = TRUE)
    r_idx = idx[1]
    c_idx = idx[2]
    print(x[idx])
    gene_name = x[r_idx, 1]
    tissues = colnames(x)[c_idx]
    print(gene_name)
    print(tissues)
    return(c(gene_name, tissues, x[idx]))
  })
  return(all_index)
}
pear_max = get_max_cor(pear_r)
pear_max_5 = get_max_cor(pear_r_5)
pear_max_10 = get_max_cor(pear_r_10)

exp1 = all_pairs.exp$CD4positivealphabetaTcell_endodermalcell[all_pairs.exp$gene_id == "TASOR2"]
his1 = lapply(all_pairs.his, function(x) x$CD4positivealphabetaTcell_endodermalcell[x$gene_id == "TASOR2"])
temp1 = as.data.frame(cbind(exp1, his1[[2]]))
temp1
# cor(temp1)
colnames(temp1) = c("exp", "his")
make_cor_plot(temp1)

prepare_cor_array <- function(gene, tissue, H_no){
  exp = all_pairs.exp[all_pairs.exp$gene_id == gene, tissue]
  tissue = gsub("trophoblastcell", "trophoblast", tissue)
  his = all_pairs.his[[H_no]]
  his = his[his$gene_id == gene, tissue]
  exp_his = as.data.frame(cbind(exp, his))
  return(exp_his)
}
make_cor_plot <- function(df, stat = TRUE){
  if (stat == TRUE){
    plot = ggplot(df, aes(x=exp, y=his)) + 
      geom_point() +
      geom_smooth(method='lm', formula= y~x) +
      stat_cor(method = "pearson", label.x = 0, label.y = -1.5, digits=2) +
      theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                         panel.background = element_blank(), aspect.ratio=1) +
      xlab(NULL) + ylab(NULL) 
  }
  else {
    plot = ggplot(df, aes(x=exp, y=his)) + 
      geom_point() +
      geom_smooth(method='lm', formula= y~x) +
      theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                         panel.background = element_blank(), aspect.ratio=1) +
      xlab(NULL) + ylab(NULL)
  }
  return(plot)
}
make_cor_plot_example <- function(gene, tissue, stat_FALSE){
  plots = list()
  for (i in 1:6){
    cor_array = prepare_cor_array(gene, tissue, i)
    if (ncol(cor_array) == 2){
      colnames(cor_array) = c("exp", "his")
      
      if (i %in% stat_FALSE) plot = make_cor_plot(cor_array, stat = FALSE)
      else plot = make_cor_plot(cor_array, stat = TRUE)
      plots[[i]] = plot
    }
    else {
      print(i)
      plots[[i]] = NULL
    }   
  }
  return(plots)
}
ex1 = make_cor_plot_example("STIM1", "H1_psoasmuscle", c(2,4,5,6))
ex2 = make_cor_plot_example("FGFR2", "CD4positivealphabetaTcell_H1", c(2,5,6))
ex3 = make_cor_plot_example("TASOR2", "CD4positivealphabetaTcell_endodermalcell", c())
all_index


get_occurence <- function(pear_df, n_points){
  all_r_vals = list()
  for (i in 1:6){
    temp = melt(pear_df[[i]][, 2:ncol(pear_df[[i]])])
    temp$H_type = histone_type_list[i]
    all_r_vals[[i]] = temp
  }
  all_r_vals = do.call(rbind, all_r_vals)
  all_r_vals = all_r_vals[!is.na(all_r_vals$value), ]
  all_r_vals$points = n_points
  return(all_r_vals)
}
pear_r_occ = get_occurence(pear_r, "1")
pear_r_occ_5 = get_occurence(pear_r_5, "5")
pear_r_occ_10 = get_occurence(pear_r_10, "10")
all_pear_r_occ = rbind(pear_r_occ, pear_r_occ_5, pear_r_occ_10)
tail(all_pear_r_occ)

cumulative_r = ggplot(all_pear_r_occ %>% 
                        group_by(H_type) %>% 
                        arrange(value) %>% 
                        mutate(rn = row_number())) + 
  geom_step(aes(x=value, y=rn, color=H_type, linetype = points)) + 
  scale_linetype_manual(values=c("solid", "dotted", "twodash")) +
  xlab("Pearson correlation R values") + ylab("Occurences") +
  theme_bw() + theme(aspect.ratio=1)
# cumulative_r
cumulative_r_ecdf = ggplot(all_pear_r_occ, aes(value, colour = H_type, linetype = points)) + 
  stat_ecdf(aes(colour = H_type)) +
  xlab("Pearson correlation R values") + ylab("Probability") +
  scale_linetype_manual(values=c("solid", "dotted", "twodash")) +
  theme_bw() + theme(aspect.ratio=1)

plot1 = ggarrange(ex1[[1]], ex1[[2]], ex1[[3]], ex1[[4]], ex1[[5]], ex1[[6]], ncol=3, nrow=2,
                  labels = histone_type_list, align = "hv")
plot1 = annotate_figure(plot1,
                        bottom = text_grob("Differential exon usage (Exon dispersion statistic)"),
                        left = text_grob("Differential histone modification (M-value)", rot = 90))
plot2 = ggarrange(ex2[[1]], ex2[[2]], ex2[[3]], ex2[[4]], ex2[[5]], ex2[[6]], ncol=3, nrow=2,
                  labels = histone_type_list, align = "hv")
plot2 = annotate_figure(plot2,
                        bottom = text_grob("Differential exon usage (Exon dispersion statistic)"),
                        left = text_grob("Differential histone modification (M-value)", rot = 90))
plot3 = ggarrange(ex3[[1]], ex3[[2]], ex3[[3]], ex3[[4]], ex3[[5]], ex3[[6]], ncol=3, nrow=2,
                  labels = histone_type_list, align = "hv")
plot3 = annotate_figure(plot3,
                        bottom = text_grob("Differential exon usage (Exon dispersion statistic)"),
                        left = text_grob("Differential histone modification (M-value)", rot = 90))
all_plots1 = ggarrange(plot1, plot2, plot3, cumulative_r,
                      nrow=2, ncol=2, align = "hv")
all_plots2 = ggarrange(plot1, plot2, plot3, cumulative_r_ecdf,
                       # labels = c("STIM1", "FGFR2", "RPS18"), 
                       nrow=2, ncol=2, align = "hv")
tiff("new3.tiff", width = 15, height = 12, units = "in", res = 200) #save pdf 20*8
all_plots1
dev.off()

tiff("new4.tiff", width = 15, height = 12, units = "in", res = 200) #save pdf 20*8
all_plots2
dev.off()


tiff("new.tiff", width = 4, height = 4, units = "in", res = 100) #save pdf 20*8
cumulative_r_ecdf
dev.off()

