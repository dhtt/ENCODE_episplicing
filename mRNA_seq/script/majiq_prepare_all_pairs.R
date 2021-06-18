library(data.table)
library(dplyr)
bam_folder = '/home/dhthutrang/ENCODE/mRNA_seq/raw_data/1'
setwd("/home/dhthutrang/ENCODE/mRNA_seq/script")
mrna_metadata = fread("/home/dhthutrang/ENCODE/mRNA_seq/raw_data/metadata.tsv")

#submeta = mrna_metadata[, c("V1", "V10")]
#colnames(submeta) = c("ID", "tissue")


submeta = mrna_metadata %>%
  mutate(tissue = gsub('[[:punct:]]|\ ', '', `Biosample term name`),
	 filename = paste(bam_folder, '/', tissue, '_', `File accession`, '.sam', sep =''),
	 rename = paste(tissue, '_', `File accession`, '.bam', sep='')
         ) %>%
  dplyr::select(tissue, `File accession`) %>%
  dplyr::group_by(tissue) %>%
  dplyr::mutate(file_acc = paste(tissue, `File accession`, sep='_'),
		group_build = paste(tissue, '=', paste(file_acc, collapse=','), sep=''),
		group_deltapsi = paste(sapply(file_acc, function(x) paste(x, '.majiq', sep='')), collapse=' ')
		) %>%
  dplyr::select(group_deltapsi) %>%
  unique()

#CHANGE NAMES BAM FILES TO TISSUE_BAM
#submeta = mrna_metadata %>%
#  mutate(tissue = gsub('[[:punct:]]|\ ', '', `Biosample term name`),
#	 command = paste('mv ', paste(bam_folder, '/', `File accession`, '.bam', sep=''), ' ', 
#			 paste(bam_folder, '/', tissue, '_', `File accession`, '.bam', sep =''), sep='')
#         ) 
#print(submeta)
#fwrite(submeta['command'], "/home/dhthutrang/ENCODE/mRNA_seq/script/name_BAM_MAJIQ.sh", row.names = FALSE, col.names = FALSE, quote = FALSE)


#GET UNIQUE TISSUES
tissue = mrna_metadata$`Biosample term name`
tissue = sapply(tissue, function(x) gsub('[[:punct:]]|\ ', '', x))
tissue = tissue[order(tissue)]
tissue = unique(tissue)
# fwrite(tissue, "/home/dhthutrang/ENCODE/mRNA_seq/script/epi_ids.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
#fwrite(tissue, "/Users/dhthutrang/Documents/BIOINFO/Episplicing/ENCODE_episplicing/mRNA_seq/script/epi_ids.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)




#group_list = submeta$group_build
group_list = submeta$group_deltapsi
for (i in tissue){
        for (j in tissue){
                if (i != j & i < j){
			name_i = group_list[grep(i, group_list)]
			name_j = group_list[grep(j, group_list)]
			print(paste('cd ', i, '_', j, ' && -grp1 ', name_i, ' -grp2 ', name_j, ' && cd .. && ls', sep=''))
                 	#lapply(c("[info]",
			#	 "bamdirs=/home/dhthutrang/ENCODE/mRNA_seq/raw_data/1",
			#	 "genome=hg38",
			#	 "[experiments]",
			#	 name_i, 
			#	 name_j), 
			#	write, 
			#	paste("majiq_config/", i, "_", j, ".txt", sep=''), 
			#	append=TRUE, ncolumns=1000)

			#python3 $MAJIQ build -j 8 -o tissue1_tissue2 transcripts >>> -c tissue1_tissue2.txt
			#write(paste("python3 $MAJIQ build /home/dhthutrang/ENCODE/refgen/majiq/hg38.ncbiRefSeq.gtf -j 8 -o ", 
			#      paste("majiq_build_res/", i, "_", j, sep=""),
			#      " -c majiq_config/", i, "_", j, ".txt", sep=""), 
			#      "majiq_build.sh",
			#      append=TRUE
			#)

			write('',
			     "majiq_deltapsi.sh",
			     append=TRUE)
                
                }
        }
}


