library(data.table)
library(dplyr)

setwd("/Users/trangdo/Documents/Episplicing/ENCODE_episplicing")
peak_metadata = fread("/Users/trangdo/Documents/Episplicing/ENCODE_episplicing/chip_seq/peak_files/metadata.tsv")
# peak_metadata[peak_metadata$V10 == "spleen" & peak_metadata$V29 == "2016-11-03", "V10"] = "spleen2016"

# submeta = peak_metadata[, c("V1", "V7", "V10", "V22")]
# colnames(submeta) = c("peak", "read", "tissue", "histone")

submeta = peak_metadata %>%
	  group_by(`Biosample term name`, `Experiment target`) %>%
	    mutate(histone_type = strsplit(`Experiment target`, '-')[[1]][1],
		            tissue_type = gsub('[[:punct:]]|\ ', '', `Biosample term name`),
			             file_names = paste(histone_type, '_', tissue_type, '.bed', sep=''), 
			             peak = paste('/home/dhthutrang/ENCODE/chip_seq/peak_files/merged/', file_names, sep=''),
				              read = paste('/home/dhthutrang/ENCODE/chip_seq/alignment_files/bed/', file_names, sep=''),
				              read = gsub('trophoblast', 'trophoblastcell', read),
					               tissue_type = gsub('trophoblast', 'trophoblastcell', tissue_type),
					               combined_files = paste(peak, read, sep="\t")) %>%
  ungroup() %>%
    dplyr::select(peak, read, tissue_type, histone_type) %>%
      unique()

for (h_type in unique(submeta$histone_type)){
	  temp = submeta %>%
		      filter(histone_type == h_type)
	        print(h_type)
		  print(paste('profile_bins --peaks=', paste(temp$peak, collapse = ','), 
			              ' --reads=', paste(temp$read, collapse = ','),
				              ' --labs=', paste(temp$tissue_type, collapse = ','),
				              ' -n multiple_testing_', paste(h_type),
					              sep=''))
}
