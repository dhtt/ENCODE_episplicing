start_time <- Sys.time()
#===== LOAD PACKAGES ======
library('stringr', quietly=TRUE)
library("data.table", quietly=TRUE)
library("dplyr", quietly=TRUE)
library("DEXSeq", quietly=TRUE)
library("optparse", quietly=TRUE)
library("BiocParallel")
suppressPackageStartupMessages( library( "DEXSeq" ) )

option_list = list(
  make_option(c("-f", "--countfolder"), type="character", default="/home/dhthutrang/ENCODE/mRNA_seq/dexseqcount/correction/count",
              help="path to folder of counts", metavar="character"),
  make_option(c("-g", "--referencegenome"), type="character", default="/home/dhthutrang/ENCODE/refgen/reference_genome.gtf",
              help="path to flattened reference genome", metavar="character"),
  make_option(c("-n", "--numcores"), type="integer", default=1,
              help="number of processing cores", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

#===== PREPARE DATA =====
setwd(opt$countfolder)
inDir = normalizePath(getwd())
count_files = list.files(inDir, pattern='*count.txt$', full.names=TRUE)
print("====================================")
print("count_files")
print(count_files)
file_names = as.data.table(str_split_fixed(basename(count_files), "\\_", 3))

gtf_files = opt$referencegenome
# log_name = file(paste(paste(epi_id1, epi_id2, sep='_'), "log", sep='.'), open = "wt")
cores = MulticoreParam(opt$numcores)

# sink(log_name, type = c("output", "message"))


# cat(paste("---> Working folder: ", opt$countfolder, sep=''), append = TRUE)
# cat("\n---> Count files: ", append = TRUE)
# cat(basename(count_files), append = TRUE)
# cat(paste("\n---> Reference genome: ", gtf_files, sep=''), append = TRUE)

sampleTable = data.frame(
  row.names = c(file_names$V2),
  condition = c(file_names$V1))
print("====================================")
print("sampleTable")
print(sampleTable)
#===== RUN DEXSEQ =====
# cat("\n---> Inputting to DEXSeq", append = TRUE)
#dxd = DEXSeqDataSetFromHTSeq(
#  count_files,
#  sampleData = sampleTable,
#  design = ~ sample + exon + condition:exon,
#  flattenedfile= normalizePath(gtf_files)
#)
#saveRDS(dxd, '/home/dhthutrang/ENCODE/mRNA_seq/dexseqcount/correction/dxd_90.RDS')
dxd = readRDS('/home/dhthutrang/ENCODE/mRNA_seq/dexseqcount/correction/dxd_90.RDS')
print(paste('dxd done: ', Sys.time()))

# cat("\n---> Getting DEXSeq result", append = TRUE)
dxd.res = DEXSeq(dxd, quiet = FALSE, BPPARAM=cores)
print(paste("dxd.res done: ", Sys.time()))
#save(dxd.res, file="/home/dhthutrang/ENCODE/mRNA_seq/dexseqcount/correction/dxd.res_90.RDS")
saveRDS(dxd.res, '/home/dhthutrang/ENCODE/mRNA_seq/dexseqcount/correction/res_90perc/dxd.res_90.RDS')
# dxd.res = readRDS('/home/dhthutrang/ENCODE/mRNA_seq/dexseqcount/correction/dxd.res.RDS')
print(paste('dxd.res saved: ', Sys.time()))

#===== SAVING RESULTS =====
# cat("\n---> Saving DEXSeq normalized counts", append = TRUE)
# dxd.count = data.frame(cbind(dxd.res[c(1,2)], counts(dxd.res, normalized = TRUE)))
# colnames(dxd.count) = c("groupID", "featureID", paste(file_names$V1, file_names$V2, sep='_'))
# normedcount_name = paste(paste(epi_id1, epi_id2, sep='_'), "normedcount.csv", sep='_')
# write.table(dxd.count, normedcount_name, quote=FALSE, sep="\t", dec=".", row.names=FALSE, col.names=TRUE)
# # dxd.count = read.csv("temp_count.csv", header=TRUE, sep = ",")
# 
# # cat("\n---> Saving DEXSeq result", append = TRUE)
# r_data_name = paste(paste(epi_id1, epi_id2, sep='_'), "RData", sep='.')
# save(dxd.res, file=r_data_name)
# 
# result_name = paste(paste(epi_id1, epi_id2, sep='_'), "res.csv", sep='_')
# write.table(as.data.frame(dxd.res[c(1,2,3,5,6,7,10)]), result_name,
#             quote=FALSE, sep="\t", dec=".", row.names=FALSE, col.names=TRUE)
# # dxd.res = read.csv(result_name, header=TRUE, sep = ",")
# 
# #print("---> Exporting HTML DEXSeq result")
# #html_name = paste(paste(epi_id1, epi_id2, sep='_'), "html", sep='_')
# #DEXSeqHTML(dxd.res, 
# #           path = html_name,
# #           FDR=0.05, color=c("#FF000080", "#0000FF80"),
# #           BPPARAM = cores)
# 
# # cat("\n===> FINISHED!", append = TRUE)
# # end_time <- Sys.time()
# # cat(paste("\nTotal time:", end_time - start_time, sep = ' '), append = TRUE)
# 
# 
# # sink()

