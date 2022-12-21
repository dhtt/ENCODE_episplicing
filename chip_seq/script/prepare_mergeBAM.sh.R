library(data.table)
library(dplyr)


# ==== Preparation ==== 
# Parse arguments for path to metadata/submetadata containing sample ID and experiment ID for grouping  
args = commandArgs(trailingOnly=TRUE)
if (length(args) == 0) {
  stop("Please input path to metadata/submetadata.tsv for grouping", call.=FALSE)
} 
metadata_path <- paste(args[1], 'submetadata.tsv', sep='/')
metadata <- fread(metadata_path, header=TRUE)


# ==== Merging BAM alignment files ==== 
# Merge by sample name and histone type
merge_files_list <- metadata %>%
  group_by(`Biosample term name`, `Experiment target`) %>%
  mutate(
    # Histone type to group the samples
    histone_type = strsplit(`Experiment target`, "-")[[1]][1],
    
    # Tissue type to group the samples
    tissue_type = gsub('[[:punct:]]|\ ', "", `Biosample term name`),

    # Merge BAM files using SAMtools 'merge' command
    dup = paste("samtools merge -@ 8", paste("merged/", histone_type, "_", tissue_type, ".bam", sep = ''), 
          paste(paste(`File accession`, '.bam', sep=''), collapse = ' ')), sep=" ") %>%
  ungroup() %>%
  dplyr::select(dup) %>%
  unique()
fwrite(merge_files_list, paste(args[1], "mergeBAM.sh", sep='/'), col.names = FALSE)
