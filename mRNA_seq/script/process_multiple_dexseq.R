library("Cairo")

res_path = "/Users/trangdo/Documents/Episplicing/ENCODE_episplicing/mRNA_seq/res_10perc/dxd.res.RDS"
html_path = "/Users/trangdo/Documents/Episplicing/ENCODE_episplicing/mRNA_seq/res_10perc/dxd.res.html"
dxd.res = readRDS(res_path)
res_df = data.frame(dxd.res)
colnames(res_df)
all_exons = paste(res_df$groupID, res_df$featureID, sep=':')
sig_exons = unique(res_df[res_df$padj < 0.05 & !is.na(res_df$padj), c('groupID', 'featureID')])
sig_exons = paste(sig_exons$groupID, sig_exons$featureID, sep=':')
length(all_exons)
length(sig_exons)
setdiff(all_exons, sig_exons)

DEXSeqHTML(dxd.res, genes='ACD', 
           path = html_path,
           FDR=0.05)
