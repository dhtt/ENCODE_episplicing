library(ggplot2)
setwd("/Users/dhthutrang/Documents/BIOINFO/Episplicing/ENCODE_episplicing/flank")

all_pairs.exp = readRDS("all_pairs.exp.RDS")
all_pairs.his = readRDS("all_pairs.his_list.RDS")

exp1 = all_pairs.exp$H1_psoasmuscle[all_pairs.exp$gene_id == "STIM1"]
his1 = lapply(all_pairs.his, function(x) x$H1_psoasmuscle[x$gene_id == "STIM1"])

exp2 = all_pairs.exp$CD4positivealphabetaTcell_H1[all_pairs.exp$gene_id == "FGFR2"]
his2 = lapply(all_pairs.his, function(x) x$CD4positivealphabetaTcell_H1[x$gene_id == "FGFR2"])

plots1 = list()
plots2 = list()

make_cor_plot <- function(df, stat = TRUE){
  if (stat == TRUE){
    plot = ggplot(df, aes(x=exp, y=his)) + 
      geom_point() +
      geom_smooth(method='lm', formula= y~x) +
      stat_cor(method = "pearson", label.x = 0, label.y = -1.5, digits=2) +
      theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                         panel.background = element_blank()) +
      xlab(NULL) + ylab(NULL) 
  }
  else {
    plot = ggplot(df, aes(x=exp, y=his)) + 
      geom_point() +
      geom_smooth(method='lm', formula= y~x) +
      theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                         panel.background = element_blank()) +
      xlab(NULL) + ylab(NULL)
  }
  return(plot)
}
for (i in 1:6){
  temp1 = as.data.frame(cbind(exp1, his1[[i]]))
  temp2 = as.data.frame(cbind(exp2, his2[[i]]))
  colnames(temp1) = c("exp", "his")
  colnames(temp2) = c("exp", "his")
  
  
  if (i %in% c(2,4,5,6)) plot1 = make_cor_plot(temp1, stat = FALSE)
  else plot1 = make_cor_plot(temp1, stat = TRUE)
  if (i %in% c(2,5,6)) plot2 = make_cor_plot(temp2, stat = FALSE)
  else plot2 = make_cor_plot(temp2, stat = TRUE)
  
  plots1[[i]] = plot1
  plots2[[i]] = plot2
}

figure <- ggarrange(ggarrange(plots1[[1]], plots1[[2]], plots1[[3]], plots1[[4]], plots1[[5]], plots1[[6]],
                              labels = c("A", "B", "C", "D", "E", "F"), align = "hv", nrow=2, ncol=3),
                    ggarrange(plots2[[1]], plots2[[2]], plots2[[3]], plots2[[4]], plots2[[5]], plots2[[6]],
                              labels = c("A", "B", "C", "D", "E", "F"), align = "hv", nrow=2, ncol=3),
                    nrow=2)
figure = annotate_figure(figure,
                bottom = text_grob("Differential exon usage (Exon dispersion statistic)"),
                left = text_grob("Differential histone modification (M-value)", rot = 90))
tiff("example_cor.tiff", width = 9, height = 9, units = "in", res = 200) #save pdf 20*8
figure
dev.off()





