## ---------------------------
##
## Script name: prepare_all_pairs.R
##
## Purpose of script: Extract the tissue names to rename SAM files for readability and generate the tissue pairs list
##
## Author: Trang Do
##
## Copyright (c) Trang Do, 2018
## Email: dhttrang@bioinformatik.uni-saarland.de
##
## ---------------------------
##
## Notes:
## - Description: Given a metadata where the tissue names, experiment and file names are stored, the names are 
##   extracted and used to rename the SAM files for easier identification. Furthermore, the list of tissue pairs are 
##   generated and splitted into chunks of pairs. These chunks are stored in separate txt files and DEXSeq analysis is
##   performed with the chunks in parallel
##
## ---------------------------


library(data.table)
library(dplyr)
library(stringr)
library(optparse, quietly = TRUE)

# ==== PREPARATION ====
# Parse arguments for input/output settings
option_list <- list(
  make_option("--SAMfolder",
    type = "character",
    help = "path to the folder of SAM files",
    default = "mRNA_seq/raw_data/sam",
    metavar = "character"
  ),
  make_option("--metadata",
    type = "character",
    help = "path to metadata.tsv",
    default = "mRNA_seq/raw_data/metadata.tsv",
    metavar = "character"
  )
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
sam_folder = opt$SAMfolder
mrna_metadata = fread(opt$metadata)

# Set working directory and input paths
working_dir <- getwd()


# ==== RENAME SAM FILES ====
submeta = mrna_metadata %>%
  mutate(tissue = gsub('[[:punct:]]|\ ', '', `Biosample term name`),
	 filename = paste(sam_folder, '/', tissue, '_', `File accession`, '.sam', sep =''),
	 rename = paste('mv ', paste(sam_folder, '/', `File accession`, '.sam', sep=''), filename)
         ) %>%
  dplyr::select(rename)
fwrite(submeta, paste(working_dir, "mRNA_seq/script/renameSAMfiles.sh", sep = "/"), 
  row.names = FALSE, col.names = FALSE, quote = FALSE)

# ==== SPLIT TISSUES PAIRS INTO CHUNK FOR PARALLEL DEXSEQ ====
tissue = mrna_metadata$`Biosample term name`
tissue = unique(sapply(tissue, function(x) gsub('[[:punct:]]|\ ', '', x)))
tissue_pairs = data.table(t(combn(tissue, 2)))
tissue_pairs = apply(tissue_pairs, 1, function(x) paste(x, collapse = ' '))
tissue_pairs = tissue_pairs[order(tissue_pairs)]

no_chunks = 11
tissue_pairs_chunks = split(tissue_pairs, ceiling(seq_along(tissue_pairs)/no_chunks))

for (i in seq_along(tissue_pairs_chunks)){
  chunk = tissue_pairs_chunks[i]
  chunk_filename = paste(paste(working_dir, "mRNA_seq/script/test/chunk_", sep = "/"),  
                         str_pad(i, 2, pad = '0'), sep="")
  fwrite(chunk, chunk_filename, row.names = FALSE, col.names = FALSE, quote = FALSE)
}