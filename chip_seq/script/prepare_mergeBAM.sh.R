library(data.table)
library(dplyr)
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  #args[1] = "out.txt"
	print(args[1])
}
file_path = paste(args[1], 'submetadata.tsv', sep='/')
file = fread(file_path, header=FALSE)
print(head(file))
rep_list = file %>%
	  group_by(V7) %>%
	    mutate(
		       dup = paste("samtools merge -@ 4", paste("merged/", V7, ".bam", sep = ''), paste(paste(V1, '.bam', sep=''), collapse = ' ')), sep=" ") %>%
  ungroup() %>%
    dplyr::select(dup) %>%
      unique()
fwrite(rep_list, paste(args[1], "mergeBAM.sh", sep='/'), col.names = FALSE)
