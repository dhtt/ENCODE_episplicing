library("data.table")
library("dplyr")
library("optparse", quietly=TRUE)


# ==== Preparation ==== 
# Parse in Epispliced project path
option_list <- list(
  make_option(c("-pd", "--projectdir"), type="character",
              help="path to Episliced project", metavar="character")
			  )
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Set woking directory to Epispliced project
setwd(opt$projectdir)


# ==== Write MAnorm/MAnorm2 arguments to file ====
# Define paths in for chip-seq analysis
chipseq_dir <- paste(opt$projectdir, "chip_seq", sep='/')
peak_dir <- paste(chipseq_dir, "peak_files", sep='/')
alignment_dir <- paste(chipseq_dir, "alignment_files", sep='/')
peak_metadata <- fread(paste(peak_dir, "metadata.tsv", sep='/'))


# Prepare data to run MAnorm/MAnorm2 script
submeta <- peak_metadata %>%
	  group_by(`Biosample term name`, `Experiment target`) %>%
	    mutate(histone_type = strsplit(`Experiment target`, '-')[[1]][1],
		            tissue_type = gsub('[[:punct:]]|\ ', '', `Biosample term name`),
			             file_names = paste(histone_type, '_', tissue_type, '.bed', sep=''), 
			             peak = paste(peak_dir, 'merged', file_names, sep='/'),
				              read = paste(alignment_dir, 'bed', file_names, sep='/'),
				              read = gsub('trophoblast', 'trophoblastcell', read),
					               tissue_type = gsub('trophoblast', 'trophoblastcell', tissue_type),
					               combined_files = paste(peak, read, sep="\t")) %>%
  ungroup() %>%
    dplyr::select(peak, read, tissue_type, histone_type) %>%
      unique()


# Write MAnorm script with info on path to peak/alignment files, sample labels and histone type
for (histone_type in unique(submeta$histone_type)){
  submeta_df <- submeta[submeta$histone_type == histone_type,]
  file_names <- as.data.frame(t(combn(submeta_df$combined_files, 2, FUN=c)))
  sample_names <- as.data.frame(t(combn(submeta_df$tissue_type, 2, FUN=c)))
  submeta_df = cbind(file_names, sample_names)
  colnames(submeta_df) = c("file_1", "file_2", "sample_1", "sample_2")
  submeta_df = submeta_df %>%
    mutate(sample1 = dplyr::if_else(sample_1 < sample_2, true = sample_1, false = sample_2),
           sample2 = dplyr::if_else(sample_1 < sample_2, true = sample_2, false = sample_1),
           file1 = dplyr::if_else(sample_1 < sample_2, true = file_1, false = file_2),
           file2 = dplyr::if_else(sample_1 < sample_2, true = file_2, false = file_1)
    ) %>%
    dplyr::select(file1, file2, sample1, sample2)
  fwrite(submeta_df, paste("chip_seq", histone_type, "all_pairs.txt", sep='/'), sep="\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
}


# Write MAnorm2 script with info on path to peak/alignment files, sample labels and histone type
for (histone_type in unique(submeta$histone_type)){
	script <- submeta %>% filter(histone_type == histone_type)
	script = paste('profile_bins --peaks=', paste(temp$peak, collapse = ','), 
		' --reads=', paste(temp$read, collapse = ','),
		' --labs=', paste(temp$tissue_type, collapse = ','),
		' -n multiple_testing_', paste(h_type),
		sep='')
	print(script)
}

