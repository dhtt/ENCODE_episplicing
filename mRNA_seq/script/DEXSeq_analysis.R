## ---------------------------
##
## Script name: DEXSeq_analysis.R
##
## Purpose of script: Perform pairwise Differential Exon Usage analysis using DEXSeq R package
##
## Author: Trang Do
##
## Copyright (c) Trang Do, 2018
## Email: dhttrang@bioinformatik.uni-saarland.de
##
## ---------------------------
##
## Notes:
## - Description:This script was used for pairwise DEU analysis. For multiple correction, 
##   DEXSeq_analysis_combined.R was used. Here, 2 samples are compared against each other using 
##   DEXSeq workflow. The script was called from execute_dexseq.sh for all possible sample pairs
##
## ---------------------------


start_time <- Sys.time()
library("stringr", quietly = TRUE)
library("data.table", quietly = TRUE)
library("dplyr", quietly = TRUE)
library("DEXSeq", quietly = TRUE)
library("optparse", quietly = TRUE)
library("BiocParallel", quietly = TRUE)


# ==== PREPARATION ====
# Parse arguments for input/output settings
option_list <- list(
  make_option(c("-f", "--countfolder"),
    type = "character",
    help = "path to folder of counts. $ENCODE_EXP/dexseqcount was used", metavar = "character"
  ),
  make_option(c("-a", "--epigenome1"),
    type = "character", default = NULL,
    help = "ID of first epigenome", metavar = "character"
  ),
  make_option(c("-b", "--epigenome2"),
    type = "character", default = NULL,
    help = "ID of second epigenome", metavar = "character"
  ),
  make_option(c("-g", "--referencegenome"),
    type = "character",
    help = "path to flattened reference genome. $REFGEN/reference_genome.gtf was used", metavar = "character"
  ),
  make_option(c("-n", "--numcores"),
    type = "integer", default = 1,
    help = "number of processing cores", metavar = "character"
  )
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

setwd(opt$countfolder)
inDir <- normalizePath(getwd())
epi_id1 <- opt$epigenome1
epi_id2 <- opt$epigenome2
pair <- paste(paste("^", epi_id1, ".*count.txt$", sep = ""), paste("^", epi_id2, ".*count.txt$", sep = ""), sep = "|")
count_files <- list.files(inDir, pattern = pair, full.names = TRUE)
file_names <- as.data.table(str_split_fixed(basename(count_files), "\\_", 3))
gtf_files <- opt$referencegenome
cores <- MulticoreParam(opt$numcores)

# Write logs
log_name <- file(paste(paste(epi_id1, epi_id2, sep='_'), "log", sep='.'), open = "wt")
sink(log_name, type = c("output", "message"))

cat(paste("---> Working folder: ", opt$countfolder, sep=''), append = TRUE)
cat("\n---> Count files: ", append = TRUE)
cat(basename(count_files), append = TRUE)
cat(paste("\n---> Reference genome: ", gtf_files, sep=''), append = TRUE)


# ===== RUN DEXSEQ =====
# Prepare sample table
cat("\n---> Preparing sample table", append = TRUE)
sampleTable <- data.frame(
  row.names = c(file_names$V2),
  condition = c(file_names$V1)
)


# Start DEXSeq analysis
cat("\n---> Inputting to DEXSeq", append = TRUE)
dxd <- DEXSeqDataSetFromHTSeq(
  count_files,
  sampleData = sampleTable,
  design = ~ sample + exon + condition:exon,
  flattenedfile = normalizePath(gtf_files)
)

cat("\n---> Getting DEXSeq result", append = TRUE)
dxd.res <- DEXSeq(dxd, quiet = FALSE, BPPARAM = cores)


# Save DEXSeq results
cat("\n---> Saving DEXSeq normalized counts", append = TRUE)
dxd.count <- data.frame(cbind(dxd.res[c(1, 2)], counts(dxd.res, normalized = TRUE)))
colnames(dxd.count) <- c("groupID", "featureID", paste(file_names$V1, file_names$V2, sep = "_"))
normedcount_name <- paste(paste(epi_id1, epi_id2, sep = "_"), "normedcount.csv", sep = "_")
write.table(dxd.count, normedcount_name, quote=FALSE, sep="\t", dec=".", row.names=FALSE, col.names=TRUE)

cat("\n---> Saving DEXSeq result", append = TRUE)
r_data_name <- paste(
  "/home/dhthutrang/ENCODE/mRNA_seq/dexseqcount/Rdata",
  paste(paste(epi_id1, epi_id2, sep = "_"), "RData", sep = "."),
  sep = "/"
)
save(dxd.res, file = r_data_name)

result_name <- paste(paste(epi_id1, epi_id2, sep = "_"), "res.csv", sep = "_")
write.table(as.data.frame(dxd.res[c(1,2,3,5,6,7,10)]), result_name,
            quote=FALSE, sep="\t", dec=".", row.names=FALSE, col.names=TRUE)
dxd.res = read.csv(result_name, header=TRUE, sep = ",")

print("---> Exporting HTML DEXSeq result")
html_name = paste(paste(epi_id1, epi_id2, sep='_'), "html", sep='_')
DEXSeqHTML(dxd.res,
          path = html_name,
          FDR=0.05, color=c("#FF000080", "#0000FF80"),
          BPPARAM = cores)

cat("\n===> FINISHED!", append = TRUE)
end_time <- Sys.time()
cat(paste("\nTotal time:", end_time - start_time, sep = ' '), append = TRUE)

sink()
