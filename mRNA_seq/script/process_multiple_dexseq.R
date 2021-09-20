library("Cairo")
setwd("/Users/trangdo/Documents/Episplicing/ENCODE_episplicing/mRNA_seq/script")

# res_path_90 = "/Users/trangdo/Documents/Episplicing/ENCODE_episplicing/mRNA_seq/res_90perc/dxd.res_90.RDS"
# html_path_90 = "/Users/trangdo/Documents/Episplicing/ENCODE_episplicing/mRNA_seq/res_90perc/dxd.res_90.html" 
all_DEU_exons = readRDS("all_DEU_genes.RDS") 
all_DEU_exons = unique(gsub(';', ':', all_DEU_exons))
all_DEU_genes = unique(sapply(all_DEU_exons, function(x) strsplit(x, split=':')[[1]][1]))
all_DEU_genes = all_DEU_genes[grep("+", all_DEU_genes, invert = T, fixed = T)]
all_DEU_genes = all_DEU_genes[grep("_", all_DEU_genes, invert = T, fixed = T)]
no_deu_gene =  length(all_DEU_genes)
no_deu_exon =  length(all_DEU_exons)

res_path_10 = "/Users/trangdo/Documents/Episplicing/ENCODE_episplicing/mRNA_seq/res_90perc/dxd.res_90.RDS"
# html_path_10 = "/Users/trangdo/Documents/Episplicing/ENCODE_episplicing/mRNA_seq/res_10perc/dxd.res.html"
dxd.res = readRDS(res_path_10)
res_df = data.frame(dxd.res)
# all_exons = rownames(res_df)

get_all_genes = function(df){
  #' Check number of all genes used in DEXSeq global
  all_genes_ = unique(df$groupID)
  print(paste("Number of investigated genes:", length(all_genes_)))
  print(head(all_genes_, 10))
  return(all_genes_)
}
all_corrected_genes = get_all_genes(res_df)
no_corrected_genes = length(all_corrected_genes)

get_sig_exon = function(df, file_name){
  #' Get significant DEU from global correction
  sig_exon_ = unique(df[df$padj < 0.05 & !is.na(df$padj), c('groupID', 'featureID')])
  sig_gene_ = unique(sig_exon_$groupID)
  sig_exon_ = paste(sig_exon_$groupID, sig_exon_$featureID, sep=':')
  no_sig_exon = length(sig_exon_)
  no_sig_gene = length(sig_gene_)
  
  #Exon
  glob_sep_overlap = intersect(sig_exon_, all_DEU_exons)
  glob_sep_union = union(sig_exon_, all_DEU_exons)
  glob_vs_sep = setdiff(sig_exon_, all_DEU_exons)
  sep_vs_glob = setdiff(all_DEU_exons, sig_exon_)
  #Gene
  glob_sep_overlap_gene = intersect(sig_gene_, all_DEU_genes)
  glob_vs_sep_gene = setdiff(sig_gene_, all_DEU_genes)
  sep_vs_glob_gene = setdiff(all_DEU_genes, sig_gene_)
  
  print(paste("Number of significant exon:", no_sig_exon))
  print(paste("Number of significant gene:", no_sig_gene))
  print(paste("Number of corrected exon:", no_deu_exon))
  print(paste("Number of corrected gene:", no_deu_gene))
  
  print(paste("DEU rate exon (Global sig/Pooled sig):", no_sig_exon/no_deu_exon))
  print(paste("DEU rate gene (Global sig/Pooled sig):", no_sig_gene/no_deu_gene))
  print(paste("Overlap - In glob - In Pooled", 
               paste(length(glob_sep_overlap), length(glob_vs_sep), length(sep_vs_glob), sep=' - ')))
  print(paste("Overlap/Union:", length(glob_sep_overlap)/length(glob_sep_union)))
  saveRDS(sig_exon_, file_name)
  return(glob_sep_overlap)
}
sig_exons = get_sig_exon(res_df, "multiple_correction_sig_exon_75low.RDS")

DEXSeqHTML(dxd.res, genes='ACD', 
           path = html_path,
           FDR=0.05)

DEXSeqHTML(dxd.res, genes='FGFR2', 
           path = html_path,
           FDR=0.05)
