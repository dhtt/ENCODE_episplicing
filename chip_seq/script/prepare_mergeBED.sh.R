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
colnames(file)
#MERGE USING TISSUE TYPE AND HISTONE TYPE
rep_list = file %>%
  group_by(`Biosample term name`, `Experiment target`) %>%
  mutate(
    nfiles = length(unique(`File accession`)),
    histone_type = strsplit(`Experiment target`, '-')[[1]][1],
    tissue_type = gsub('[[:punct:]]|\ ', '', `Biosample term name`),
    file_name = paste('merged/', histone_type, '_', tissue_type, '.bed', sep = ''),
    dup = dplyr::if_else(nfiles == 1, 
                         true = paste('cp', paste(`File accession`, '.bed', sep=''), file_name, sep =' '),
                         false = paste('cat', 
                                       paste(paste(`File accession`, '.bed', sep=''), collapse = ' '),
                                       '>', file_name, sep=' ')),
    sort = paste('sort -k1,1 -k2,2n ', file_name ,' > ', paste(file_name, '.sorted.bed', sep='')), 
    merge = dplyr::if_else(nfiles == 1,
                           true = paste('mv ', paste(file_name, '.sorted.bed ', sep=''), file_name,sep = ''),
                           false = paste('bedtools merge -c 4,5,6,7 -o collapse,sum,distinct,sum -i ', 
                                         paste(file_name, '.sorted.bed', sep=''), ' > ', file_name,sep = '')),
    remove = dplyr::if_else(nfiles == 1,
                            true = paste('rm ', paste(file_name, '.sorted.bed', sep=''), sep =''), 
                            false = paste('')), 
    all = paste(dup, sort, merge, remove, sep =';')) %>%
  ungroup() %>%
  dplyr::select(all) %>%
  unique()
fwrite(rep_list, paste(args[1], 'mergeBED.sh', sep='/'), col.names = FALSE, quote = FALSE)
