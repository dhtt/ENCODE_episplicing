## ---------------------------
##
## Script name: get_dexseq_plot.R
##
## Purpose of script: Get DEXSeq-differential exon usage plot for any gene in any tissue pairs comparison. The
## script is used to generate figure 4A
##
## Author: Trang Do
##
## Copyright (c) Trang Do, 2018
## Email: dhttrang@bioinformatik.uni-saarland.de
##
## ---------------------------

#===== LOAD PACKAGES ======
library("DEXSeq", quietly = TRUE)
library("optparse", quietly = TRUE)


option_list <- list(
  make_option(c("-a", "--epigenome1"), 
    type = "character", 
    default = NULL,
    help = "ID of first epigenome", 
    metavar = "character"
    ),
  make_option(c("-b", "--epigenome2"), 
    type = "character", 
    default = NULL,
    help = "ID of second epigenome", 
    metavar="character"),
  make_option(c("-g", "--gene_ID"), 
    type = "character", 
    default = NULL,
    help = "chosen gene for plotting", 
    metavar = "character"
    ),
  make_option(c("-d", "--data_path"),
    type = "character", 
    default = "mRNA_seq/dexseqcount/Rdata",
    help = "path to folder where Rdata of the DEXseq results are stored", 
    metavar = "character"
  ),
  make_option(c("-o", "--general_analysis_results"),
    type = "character",
    help = "path to the folder where the general results from the analysis are stored and shared between processes",
    metavar = "character",
    default = "general_analysis_results"
  )
  make_option(c("-w", "--width"), 
    type = "integer", 
    default = 8,
    help = "set image width", 
    metavar="character"),
  make_option(c("-h", "--height"), 
    type = "integer", 
    default = 5,
    help = "set image height", 
    metavar="character"),
  make_option(c("-r", "--resolution"), 
    type = "integer", 
    default = 300,
    help = "set image resolution", 
    metavar="character")
)

opt_parser <- OptionParser(option_list = option_list, add_help_option = FALSE)
opt <- parse_args(opt_parser)
epi_id1 <- opt$epigenome1
epi_id2 <- opt$epigenome2
gene_ID <- opt$gene_ID
image_height <- opt$height
image_width <- opt$width
image_res <- opt$resolution
general_analysis_results <- opt$general_analysis_results
#===== Plot gene ===== 
print("..... Gather data \n")
RData_filename <- paste(paste(paste(epi_id1, epi_id2, sep = "_"), "*.RData", sep = "."),
                        paste(paste(epi_id2, epi_id1, sep = "_"), "*.RData", sep = "."), sep = "|")
RData_files <- list.files(opt$data_path, pattern = RData_filename, full.names = FALSE)

print("..... Plotting from: \n")
print(RData_files)

for (file in RData_files){
  load(file)
  plot_name <- paste(paste(general_analysis_results, file, sep = "/"), gene_ID, 'tiff', sep = ".")
  tiff(plot_name, units = "in", width = image_width, height = image_height, res = image_res)
  plotDEXSeq(dxd.res, gene_ID, legend = TRUE, FDR = 0.05, color = c("#EE442F", "#63ACBE"), lwd = 1.8,
            norCounts = TRUE, splicing = TRUE, expression = TRUE)
  dev.off()
}
