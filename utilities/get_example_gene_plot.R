library(Gviz)
library(GenomicRanges)
library(rtracklayer)
library(stringr)
library(data.table)
library(dplyr)
epigenomes_names
col.highlight = "#fffa70"
fontsize.title = 17 #normally 17
type.color = c("#EE442F", "#63ACBE")
margin = c(20, 20)
gen = "hg19"
background.title = "#5D2ED2"
col.title = "white"
#======================================================================================
gene_name = "FGFR2" #INPUT
epi_id1 = "CD4positivealphabetaTcell" #INPUT
epi_id2 = "H1" #INPUT
example_gene_folder = normalizePath(paste("/Users/dhthutrang/Documents/BIOINFO/Episplicing/ENCODE_episplicing/utilities",
                                          paste(epi_id1, epi_id2, sep ='_'), sep = '/')) #INPUT
epi_name1 = "CD4-positive alpha beta T cell" #INPUT
epi_name2 = "H1 cell" #INPUT
reverseStrand = TRUE #INPUT
# grouping = c(rep("control", 2), rep("treated", 3))
all_files = list.files(example_gene_folder, full.names = TRUE)
all_files
#===== EXP =====
transcript_file = all_files[grep("NCBI", all_files)]
transcript = import.gff(transcript_file)
transcript = as(transcript, "GRanges")
transcript = transcript[transcript$type == "exon",]

" ++++++++++++++NM_000141.5+XM_024447891.1+XM_024447891.1+NM_023029.2+NM_001144915.2+XM_017015924.2+XM_024447890.1+NM_001144914.1+NM_001144918.2"
# omit_transcript = c("NM_001382568.1", "NM_001277961.3", "NM_001382574.1", "NM_001382572.1",
#                     "NM_001382580.1", "NM_001382581.1", "NM_001382569.1", "NM_001382571.1",
#                     "NM_001382579.1", "NM_001277962.2", "NM_001382578.1") #STIM1
omit_transcript = c("XM_024447892.1", "NM_001144915.2", "NM_001144916.2",
                    "XM_024447890.1", "NM_022970.3", "NM_001144917.2",
                    "NM_023029.2", "NM_001144914.1",
                    "NM_001144919.2", "XM_024447887.1", "XM_024447888.1", "XM_024447889.1",
                    "NM_001144918.2", "XM_017015924.2","XM_017015925.2", "XM_024447891.1",
                    )
transcript = transcript[transcript$transcript %in% omit_transcript == FALSE,]

chr <- as.character(unique(seqnames(transcript)))[[1]]
itrack <- IdeogramTrack(genome = gen, chromosome = chr, name = gene_name,
                        fontsize.title = fontsize.title)
gtrack <- GenomeAxisTrack(labelPos = "above",
                          fontsize.title = fontsize.title)
grtrack <- GeneRegionTrack(transcript,
                           genome = gen, chromosome = chr,
                           name = "Transcripts", 
                           transcriptAnnotation = "transcript",
                           fontsize.title = fontsize.title,
                           col = "black", fill="#85c0f9",
                           cex.group = 0.5
)

transcript_sig_file = all_files[grep("res.csv", all_files)]
transcript_sig = fread(transcript_sig_file, header = TRUE)

#----- Prepare expression level -----
# expression_file = fread(all_files[grep("normedcount.csv", all_files)])
min_max_nnorm <- function(x){return ((x - min(x))/ (max(x) - min(x)))}
counts = subset(expression_file, select = -c(1,2))
counts = counts %>% 
  mutate_all(., as.numeric) %>%
  mutate_all(., min_max_nnorm)

expression = cbind(transcript_sig[, c(8, 10, 14, 11, 12)], counts)
colnames(expression) = c("chr", "feature", "strand", "start", "end", paste("count", seq(1, ncol(expression)-5), sep = "_"))
expression = as(expression,"GRanges")

eTrack <- DataTrack(expression,
                    legend = FALSE, name = "Exon usage",
                    groups = grouping,
                    type = c("a", "p", "confint"),
                    col = type.color)
# plotTracks(eTrack)

#----- Prepare highlighted transcript -----
transcript_sig = transcript_sig[transcript_sig$padj <= 0.05, c(6, 11 ,12)]
colnames(transcript_sig) = c("padj", "start", "end")

# transcript_sig = transcript_sig[transcript_sig$padj <= 0.05, c(6)]
# # temp = rbindc(4107707, 4107773)
# temp = rbind(c(3876894, 3876932), c(3988782, 3988912), c(4076756, 4076867), c(4080511, 4080626), c(4091256, 4091433),
#               c(4104112, 4104212), c(4104493, 4104728), c(4107707, 4107773), c(4107774, 4108091), c(4112512, 4114438))
# transcript_sig = cbind(transcript_sig, temp)
# colnames(transcript_sig) = c("padj", "start", "end")
# transcript_sig = transcript_sig[-c(6,8), ] #FOR STIM1

grtrack_sig <- HighlightTrack(trackList = grtrack, chromosome = chr,
                              start = transcript_sig$start - 100,
                              end = transcript_sig$end + 100,
                              col = "yellow", fill = "yellow", alpha = 1,
                              inBackground = TRUE)
# if (dim(transcript_sig)[[1]] != 0){
#   grtrack_sig <- HighlightTrack(trackList = list(grtrack, eTrack), chromosome = chr,
#                                 start = transcript_sig$start - 300,
#                                 end = transcript_sig$end + 300,
#                                 col = "yellow", fill = "yellow", alpha = 1,
#                                 inBackground = TRUE)
# } else {
#   grtrack_sig = list(grtrack, eTrack)
# }

plotTracks(list(itrack, gtrack, grtrack_sig), extend.left = 2500, extend.right = 2500,
           add53 = TRUE, add35 = TRUE, cex = 0.5, reverseStrand = FALSE)

#===== HIS =====
library(dplyr)
# his_type = c("H3K27ac", "H3K36me3", "H3K4me1")
his_type = c("H3K27ac", "H3K27me3", "H3K36me3", "H3K4me3", "H3K4me1", "H3K9me3")
# his_type = c("H3K27ac", "H3K27me3", "H3K36me3", "H3K4me3", "H3K9me3")
prepare_hisdiff_example <- function(epi_id1, epi_id2, gene, type){
  file_id_exon = paste(".*exon_", paste(type, "txt", sep='.'), sep='.*')
  file_id_pi = paste(".*pi_", paste(type, "txt", sep='.'), sep='.*')
  hisdiff.exon = fread(normalizePath(all_files[grep(file_id_exon, all_files, fixed = FALSE)]), 
                       header=FALSE, stringsAsFactors = FALSE, quote=FALSE)
  hisdiff.pi = fread(normalizePath(all_files[grep(file_id_pi, all_files, fixed = FALSE)]), 
                     header=FALSE, stringsAsFactors = FALSE, quote=FALSE)
  hisdiff_tables_list = list(hisdiff.exon, hisdiff.pi)
  
  for (i in c(1,2)){
    hisdiff_table = hisdiff_tables_list[[i]]
    hisdiff_table <- hisdiff_table %>%
      dplyr::filter(V3 == "exonic_part" | V3 == "promoter" | V3 == "intron") %>%
      dplyr::mutate(V11 = dplyr::if_else(V11 == -1, V4, V11), V12 = dplyr::if_else(V12 == -1, V5, V12),
             
             V15 = as.numeric(as.character(V15)), V16 = as.numeric(as.character(V16)),
             V15 = dplyr::if_else(is.na(V15), 0, V15), V16 = dplyr::if_else(is.na(V16), 0, V16),
             
             V13 = as.numeric(as.character(V13)), V14 = as.numeric(as.character(V14)),
             V13 = dplyr::if_else(is.na(V13), 0, V13), V14 = dplyr::if_else(is.na(V14), 1, V14),
             V17 = V9) %>%
      dplyr::select(c(1,3,7,11,12,13,14,15,16,17))

    colnames(hisdiff_table) = c("chr", "feature", "strand", "start", "end","m-val","p_val", epi_id1, epi_id2, "gene")
    hisdiff_tables_list[[i]] = hisdiff_table
  }
  hisdiff_tables_list = rbind(hisdiff_tables_list[[1]], hisdiff_tables_list[[2]])
  hisdiff_tables_list = as(hisdiff_tables_list, "GRanges")
  return(hisdiff_tables_list)
}
get_all_his_tables <- function(m_thres = 1){
  all_his_counts = vector("list", length(his_type))
  all_his_sigs = vector("list", length(his_type))
  for (i in 1:length(his_type)){
    type = his_type[i]
    print(type)
    his_table = prepare_hisdiff_example(epi_id1, epi_id2, gene_name, type)
    
    print("his_table.count")
    his_table.count = his_table[, c(epi_id1, epi_id2)]   
    print(his_table.count)
    print("his_table.sig")
    his_table.sig = as.data.frame((his_table[ his_table$`p_val` <= 0.05 &
                                               abs(his_table$`m-val`) >= m_thres 
                                             ]))
    print(head(his_table.sig))
    all_his_counts[[i]] = his_table.count
    all_his_sigs[[i]] = his_table.sig
  }
  # return(his_table)
  return(list(all_his_counts, all_his_sigs))
} 
all_his_counts = get_all_his_tables(m_thres = 1)
all_his_sigs = all_his_counts[[2]]
all_his_counts = all_his_counts[[1]]
lapply(all_his_sigs, function(x) dim(x)[[1]])

get_all_his_track.sigs <- function(expansion){
  all_his_track.sigs = vector("list", length(all_his_counts))
  for (i in 1:length(all_his_counts)){
    his_count = all_his_counts[[i]]
    his_sig = all_his_sigs[[i]]
    if (dim(his_sig)[[1]] != 0){
      if (i == length(all_his_counts)) {
        hTrack <- DataTrack(his_count, name = his_type[i],
                            type = c("a"), legend = TRUE,
                            groups = c(epi_name1, epi_name2),
                            col = type.color)
      } else {
        hTrack <- DataTrack(his_count, name = his_type[i],
                            type = c("a"), legend = FALSE,
                            groups = c(epi_name1, epi_name2),
                            col = type.color)
      }
      hTrack.sig <- HighlightTrack(trackList = hTrack, chromosome = chr,
                                   start = his_sig$start  - expansion,
                                   end = his_sig$end  + expansion,
                                   col = col.highlight, fill = col.highlight)
      all_his_track.sigs[[i]] = hTrack.sig
    }
    else {
      hTrack <- DataTrack(his_count, name = his_type[i],
                          type = c("a"), legend = FALSE,
                          col = type.color)
      all_his_track.sigs[[i]] = hTrack
    }
  }
  return(all_his_track.sigs)
}
all_his_track.sigs = get_all_his_track.sigs(0)


# temp = all_his_sigs[[6]]
# sum(temp$end - temp$start)/(max(end(transcript)) - min(start(transcript)))*100

#===== ALL TRACKS =====
all_tracks = list(itrack, gtrack, grtrack_sig)
all_tracks = append(all_tracks, all_his_track.sigs)
# all_tracks = append(all_tracks, mTrack_sig)
plotTracks(all_tracks,
           extend.left = 2500, extend.right = 2500,
           background.title = background.title,
           fontsize.title = fontsize.title,
           fontsize = fontsize.title,
           col.title = col.title,
           margin = margin,
           reverseStrand = reverseStrand,
           stackHeight = 0.5,
           add53 = TRUE, add35 = TRUE, cex = 0.5, lwd=1.5, cex.sampleNames = 1)


tiff(paste(paste(epi_id1, epi_id2, gene_name, sep='_'), "tiff", sep='.'), 
     units="in", width=12, height=15, res=300)
plotTracks(all_tracks,
           extend.left = 5000, extend.right = 500,
           background.title = background.title,
           fontsize.title = fontsize.title,
           # fontsize = fontsize.title,
           col.title = col.title,
           margin = margin, 
           reverseStrand = reverseStrand,
           stackHeight = 0.75, cex = 0.6, lwd=1.5, frame=TRUE, cex.axis= 0.3, 
           lineheight = 0.25, cex.legend= 0.75, fontsize.legend = fontsize.title)
dev.off()


# col.highlight = "#fffa70"
# fontsize.title = 22 #normally 17
# type.color = c("#EE442F", "#63ACBE")
# margin = c(20, 20)
# gen = "hg19"
# background.title = "#5D2ED2"
# col.title = "white"
# 
# tiff(paste(paste(epi_id1, epi_id2, gene_name, sep='_'), "tiff", sep='.'), 
#      units="in", width=8, height=12, res=300)
# plotTracks(all_tracks,
#            extend.left = 1000, extend.right = 500,
#            background.title = background.title,
#            fontsize.title = fontsize.title,
#            fontsize = fontsize.title,
#            col.title = col.title,
#            groups = c(epi_name1, epi_name2),
#            margin = margin, 
#            reverseStrand = reverseStrand,
#            stackHeight = 0.75, cex = 0.5, lwd=1.5, frame=TRUE, cex.axis= 0.3, lineheight=0.1)
# dev.off()


#===== MET =====
prepare_metdiff_example <- function(epi_id1, epi_id2, gene){
  metdiff.exon = fread(all_files[grep('.*exon_.*normedratio.csv.txt', all_files, fixed = FALSE)], 
                       header=FALSE, stringsAsFactors = FALSE, quote=FALSE)
  metdiff.pi = fread(all_files[grep('.*pi_.*normedratio.csv.txt', all_files, fixed = FALSE)], 
                     header=FALSE, stringsAsFactors = FALSE, quote=FALSE)
  head(metdiff.exon)
  head(metdiff.pi)
  metdiff_tables_list = list(metdiff.exon, metdiff.pi)
  
  for (i in range(1:2)){
    metdiff_table = metdiff_tables_list[[i]]
    metdiff_table <- metdiff_table %>%
      filter(V3 == "exonic_part" | V3 == "promoter" | V3 == "intron") %>%
      mutate(V11 = if_else(V11 == -1, V4, V11), V12 = if_else(V12 == -1, V5, V12),
             V14 = as.numeric(as.character(V14)), V15 = as.numeric(as.character(V15)),
             V14 = if_else(is.na(V14), 0, V14), V15 = if_else(is.na(V15), 0, V15),
             V16 = gene) %>%
      dplyr::select(c(1,3,7,11,12,14,15,16))
    
    colnames(metdiff_table) = c("chr", "feature", "strand", "start", "end", epi_id1, epi_id2, "gene")
    metdiff_tables_list[[i]] = metdiff_table
  }
  metdiff_tables_list = rbind(metdiff_tables_list[[1]], metdiff_tables_list[[2]])
  metdiff_tables_list = as(metdiff_tables_list, "GRanges")
  return(metdiff_tables_list)
}
met_count = prepare_metdiff_example(epi_id1, epi_id2, gene_name)
mTrack <- DataTrack(met_count, name = "Methylation",
                    groups = c(epi_name1, epi_name2),
                    type = c("a"),
                    col = type.color,
                    legend = TRUE)

get_sig_met <- function(epi_id1, epi_id2, gene, thres_metdiff){
  metdiff = fread(normalizePath(all_files[grep('.*diff.txt.txt', all_files, fixed = FALSE)]), 
                  header=FALSE, stringsAsFactors = FALSE, quote=FALSE)
  metdiff_table <- metdiff %>%
    mutate(V11 = if_else(V11 == -1, V4, V11), V12 = if_else(V12 == -1, V5, V12),
           V15 = as.numeric(as.character(V15)), V16 = as.numeric(as.character(V16)),
           V15 = if_else(is.na(V15), 1, V15), V16 = if_else(is.na(V16), 0, V16),
           V17 = gene) %>%
    dplyr::select(c(1,3,7,11,12,15,16,17))
  colnames(metdiff_table) = c("chr", "feature", "strand", "start", "end", "q-val", "metdiff","gene")
  metdiff_table = metdiff_table[metdiff_table$`q-val` <= 0.05 & abs(metdiff_table$`metdiff`) >= thres_metdiff, ]
  metdiff_table = as(metdiff_table, "GRanges")
  metdiff_table = as.data.frame(metdiff_table)
  return(metdiff_table)
}
met_sig = get_sig_met(epi_id1, epi_id2, gene_name, 25)

if (dim(met_sig)[[1]] != 0){
  mTrack_sig <- HighlightTrack(trackList = mTrack, chromosome = chr,
                               start = met_sig$start - 100,
                               end = met_sig$end + 200,
                               groups = c(epi_id1, epi_id2),
                               col = "#fffa70", fill = "#fffa70")
} else {
  mTrack_sig = mTrack
}

plotTracks(list(itrack, gtrack, grtrack_sig, mTrack_sig),
           extend.left = 2500, extend.right = 2500,
           background.title = "#597ca8")

