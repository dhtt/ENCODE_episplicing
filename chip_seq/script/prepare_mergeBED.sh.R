library(data.table)
library(dplyr)
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop('At least one argument must be supplied (input file).n', call.=FALSE)
} else if (length(args)==1) {
  # default output file
  #args[1] = 'out.txt'
  print(args[1])
}
file_path = paste(args[1], 'metadata.tsv', sep='/')
file = fread(file_path, header=TRUE)
# file = fread('/Users/dhthutrang/Documents/BIOINFO/Episplicing/ENCODE_episplicing/chip_seq/peak_files/metadata.tsv', header=TRUE)

head(file)
#MERGE USING TISSUE TYPE AND HISTONE TYPE
rep_list = file %>%
  group_by(`Biosample term name`, `Experiment target`) %>%
  mutate(
    nfiles = length(unique(`File accession`)),
    histone_type = strsplit(`Experiment target`, '-')[[1]][1],
    tissue_type = gsub('[[:punct:]]|\ ', '', `Biosample term name`),
    dup = dplyr::if_else(nfiles == 1, 
                         true = paste('mv ', paste(`File accession`, '.bed', sep=''), paste('merged/', histone_type, '_', tissue_type, '.bed', sep = '') ),
                         false = paste('bedtools merge -i', 
                                       paste(paste(`File accession`, '.bed', sep=''), collapse = ' '),
                                       '-o sum >', 
                                       paste('merged/', histone_type, '_', tissue_type, '.bed', sep = ''), sep=' '))) %>%
  ungroup() %>%
  dplyr::select(dup) %>%
  unique()
fwrite(rep_list, paste(args[1], 'mergeBED.sh', sep='/'), col.names = FALSE)
