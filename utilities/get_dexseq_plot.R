## ---------------------------
##
## Script name: get_dexseq_plot.R
##
## Purpose of script: Get DEXSeq-differential exon usage plot for any gene in any tissue pairs comparison
##
## Author: Trang Do
##
## Copyright (c) Trang Do, 2018
## Email: dhttrang@bioinformatik.uni-saarland.dels

##
## ---------------------------
##
## Notes:
## - Description: -
## - Preceeding script: -
## - Succeeding script: -
##
## ---------------------------

#===== LOAD PACKAGES ======
library("DEXSeq", quietly = TRUE)
library("optparse", quietly = TRUE)

setwd("/home/dhthutrang/ENCODE/mRNA_seq/dexseqcount/Rdata")

option_list <- list(
  make_option(c("-a", "--epigenome1"), type = "character", default = NULL,
              help = "ID of first epigenome", metavar="character"),
  make_option(c("-b", "--epigenome2"), type = "character", default = NULL,
              help = "ID of second epigenome", metavar="character"),
  make_option(c("-g", "--gene_ID"), type = "character", default = NULL,
              help = "chosen gene for plotting", metavar="character"),
  make_option(c("-w", "--width"), type = "integer", default = 8,
              help = "set image width", metavar="character"),
  make_option(c("-h", "--height"), type = "integer", default = 5,
              help = "set image height", metavar="character"),
  make_option(c("-r", "--resolution"), type = "integer", default = 300,
              help = "set image resolution", metavar="character")
)

opt_parser <- OptionParser(option_list = option_list, add_help_option = FALSE)
opt <- parse_args(opt_parser)
inDir <- normalizePath(getwd())
epi_id1 <- opt$epigenome1
epi_id2 <- opt$epigenome2
gene_ID <- opt$gene_ID
image_height <- opt$height
image_width <- opt$width
image_res <- opt$resolution


#===== Plot gene ===== 
print("..... Gather data \n")
RData_filename <- paste(paste(paste(epi_id1, epi_id2, sep = "_"), "*.RData", sep = "."),
                        paste(paste(epi_id2, epi_id1, sep = "_"), "*.RData", sep = "."), sep = "|")
RData_files <- list.files(inDir, pattern = RData_filename, full.names = FALSE)

print("..... Plotting from: \n")
print(RData_files)

for (file in RData_files){
  load(file)
  plot_name <- paste(paste("/home/dhthutrang/ENCODE/utilities/plot", file, sep = "/"), gene_ID, 'tiff', sep = ".")
  tiff(plot_name, units = "in", width = image_width, height = image_height, res = image_res)
  plotDEXSeq(dxd.res, gene_ID, legend = TRUE, FDR = 0.05, color = c("#EE442F", "#63ACBE"), lwd = 1.8,
            norCounts = TRUE, splicing = TRUE, expression = TRUE)
  dev.off()
}
