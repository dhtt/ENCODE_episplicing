## ---------------------------
##
## Script name: DEXSeq_analysis_combined.R
##
## Purpose of script: Perform combined Differential Exon Usage analysis using DEXSeq R package
##
## Author: Trang Do
##
## Copyright (c) Trang Do, 2018
## Email: dhttrang@bioinformatik.uni-saarland.de
##
## ---------------------------
##
## Notes: 
## - Description: This script was used for multiple correction of DEU analysis. For pairwise 
## comparison, DEXSeq_analysis.R was used
## - Preceeding script: - multiple_correction.R
## - Succeeding script: - process_multiple_dexseq.R
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
  make_option(c("-f", "--count_folder"),
    type = "character", 
    default = "mRNA_seq/dexseqcount",
    help = "path to folder of counts", metavar = "character"
  ),
  make_option(c("-g", "--reference_genome"),
    type = "character", 
    default = "refgen/reference_genome.gtf",
    help = "path to flattened reference genome", 
    metavar = "character"
  ),
  make_option(c("-n", "--num_cores"),
    type = "integer", 
    default = 1,
    help = "number of processing cores", 
    metavar = "character"
  )
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

count_files <- list.files(paste(opt$count_folder, "correction/count", sep = "/"), pattern = "*count.txt$", full.names = TRUE)
file_names <- as.data.table(str_split_fixed(basename(count_files), "\\_", 3))
gtf_files <- opt$reference_genome
cores <- MulticoreParam(opt$num_cores)
correction_result_path <- paste(opt$count_folder, "correction/count/res", sep = "/")

# Write logs
log_name <- file(paste(paste(epi_id1, epi_id2, sep='_'), "log", sep='.'), open = "wt")
sink(log_name, type = c("output", "message"))

cat(paste("---> Working folder: ", opt$count_folder, sep=''), append = TRUE)
cat("\n---> Count files: ", append = TRUE)
cat(basename(count_files), append = TRUE)
cat(paste("\n---> Reference genome: ", gtf_files, sep=''), append = TRUE)


# ===== RUN DEXSEQ =====
# Prepare sample table
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
saveRDS(dxd, paste(correction_result_path, "dxd.RDS", sep = "/"))
print(paste("dxd done: ", Sys.time()))

cat("\n---> Getting DEXSeq result", append = TRUE)
dxd.res <- DEXSeq(dxd, quiet = FALSE, BPPARAM = cores)
print(paste("dxd.res done: ", Sys.time()))


# Save DEXSeq results
saveRDS(dxd.res, paste(correction_result_path, "dxd.res.RDS", sep = "/"))
print(paste("dxd.res saved: ", Sys.time()))


cat("\n===> FINISHED!", append = TRUE)
end_time <- Sys.time()
cat(paste("\nTotal time:", end_time - start_time, sep = ' '), append = TRUE)
sink()
