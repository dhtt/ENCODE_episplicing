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
 

# ==== Merging BED peak files ==== 
# Merge by sample name and histone type
merge_files_list <- file %>%
  group_by(`Biosample term name`, `Experiment target`) %>%
  mutate(
    nfiles = length(unique(`File accession`)),

    # Histone type to group the samples
    histone_type = strsplit(`Experiment target`, '-')[[1]][1],

    # Tissue type to group the samples
    tissue_type = gsub('[[:punct:]]|\ ', '', `Biosample term name`),

    # Write file names for each group according to histone and tissue type
    file_name = paste('merged/', histone_type, '_', tissue_type, '.bed', sep = ''),

    # If there is only 1 sample, use the sample. Else, merge the samples by 'cat' command
    dup = dplyr::if_else(nfiles == 1, 
                         true = paste('cp', paste(`File accession`, '.bed', sep=''), file_name, sep =' '),
                         false = paste('cat', 
                                      paste(paste(`File accession`, '.bed', sep=''), collapse = ' '),
                                      '>', file_name, sep=' ')
                        ),

    # Sort BED files using BEDtools 'sort'
    sort = paste('sort -k1,1 -k2,2n ', file_name ,' > ', paste(file_name, '.sorted.bed', sep='')), 

    # Merge overlapsed reads in BED files using BEDtools 'merge'
    merge = dplyr::if_else(nfiles == 1,
                           true = paste('mv ', paste(file_name, '.sorted.bed ', sep=''), file_name,sep = ''),
                           false = paste('bedtools merge -c 4,5,6,7 -o collapse,sum,distinct,sum -i ', 
                                         paste(file_name, '.sorted.bed', sep=''), ' > ', file_name,sep = '')),

    # Remove temporary files                                    
    remove = dplyr::if_else(nfiles == 1,
                            true = paste('rm ', paste(file_name, '.sorted.bed', sep=''), sep =''), 
                            false = paste('')), 

    # Execute all commands
    all_commands = paste(dup, sort, merge, remove, sep =';')) %>%
  ungroup() %>%
  dplyr::select(all) %>%
  unique()
fwrite(merge_files_list, paste(args[1], 'mergeBED.sh', sep='/'), col.names = FALSE, quote = FALSE)
