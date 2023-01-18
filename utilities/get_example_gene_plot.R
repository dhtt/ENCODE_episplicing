## ---------------------------
##
## Script name: analysis_pipeline.R
##
## Purpose of script: Perform the main analysis on DEU-DHM correlation data
##
## Author: Trang Do
##
## Copyright (c) Trang Do, 2018
## Email: dhttrang@bioinformatik.uni-saarland.de
##
## ---------------------------
##
## Notes:
## - Description: This script plot the DEXSeq results (normalized exon usage and DEU exons), MAnorm results (normalized 
##   ChIP-seq signals and DHM exons) for an example gene in a comparison between two epigenomes. The generated plot is 
##   included in Figure 4. To gather the files needed for plotting, prepare_example_gene.sh is used.
## - Preceeding script: prepare_example_gene.sh
## - Succeeding script: -
##
## ---------------------------

#===== LOAD PACKAGES ======
library("Gviz", quietly = TRUE)
library("GenomicRanges", quietly = TRUE)
library("rtracklayer", quietly = TRUE)
library("stringr", quietly = TRUE)
library("data.table", quietly = TRUE)
library("dplyr", quietly = TRUE)

#==== PREPARATION ====
# Parse arguments for input/output settings
col.highlight <- "#fffa70"
fontsize.title <- 17
type.color <- c("#EE442F", "#63ACBE")
margin <- c(20, 20)
background.title <- "#5D2ED2"
col.title <- "white"
gene_name <- "SEPTIN9" # INPUT
epi_id1 <- "neuronalstemcell" # INPUT
epi_id2 <- "spleen" # INPUT
example_gene_folder <- normalizePath(paste("/home/dhthutrang/ENCODE/utilities",
  paste("example", epi_id1, epi_id2, gene_name, sep = "_"),
  sep = "/"
))
example_gene_folder <- normalizePath(paste("/Users/dhthutrang/Documents/BIOINFO/Episplicing/ENCODE_episplicing/utilities",
  paste(epi_id1, epi_id2, sep = "_"),
  sep = "/"
)) # INPUT
epi_name1 <- "Neuronal stem cell" # INPUT
epi_name2 <- "Spleen" # INPUT
reverseStrand <- TRUE # INPUT
all_files <- list.files(example_gene_folder, full.names = TRUE)
# omit_transcript = c("NM_001382568.1", "NM_001277961.3", "NM_001382574.1", "NM_001382572.1",
#                     "NM_001382580.1", "NM_001382581.1", "NM_001382569.1", "NM_001382571.1",
#                     "NM_001382579.1", "NM_001277962.2", "NM_001382578.1") #STIM1
omit_transcript <- c(
  "XM_024447892.1", "NM_001144915.2", "NM_001144916.2",
  "XM_024447890.1", "NM_022970.3", "NM_001144917.2", "NM_023029.2", "NM_001144914.1",
  "NM_001144919.2", "XM_024447887.1", "XM_024447888.1", "XM_024447889.1",
  "NM_001144918.2", "XM_017015924.2", "XM_017015925.2", "XM_024447891.1",
)
histone_types <- c("H3K27ac", "H3K27me3", "H3K36me3", "H3K4me3", "H3K9me3")


#==== PLOT DIFFERENTIAL EXON USAGE =====

#---- 1. Prepare transcripts and the genomic region to plot on ----

# Read reference genome file for gene's region and transcripts 
transcript_file <- all_files[grep("NCBI", all_files)]
transcript <- as(import.gff(transcript_file), "GRanges")

# Only show exons and main transcripts
transcript <- transcript[transcript$type == "exon", ]
transcript <- transcript[transcript$transcript %in% omit_transcript == FALSE, ]

#---- 2. Prepare highlighted transcript ----
# Read DEXSeq result file for exons with DEUs (FDR p-value <= 0.05)
DEU_file <- all_files[grep("res.csv", all_files)]
DEU_table <- fread(DEU_file, header = TRUE)
DEU_table <- DEU_table[DEU_table$padj <= 0.05, c(6, 11, 12)]
colnames(DEU_table) <- c("padj", "start", "end")

#---- 3. Plot genomic region with highlighted DEU exons ----
chr <- as.character(unique(seqnames(transcript)))[[1]]
genome <- "hg38"
itrack <- IdeogramTrack(genome = genome, chromosome = chr, name = gene_name, fontsize.title = fontsize.title)
gtrack <- GenomeAxisTrack(labelPos = "above", fontsize.title = fontsize.title)
grtrack <- GeneRegionTrack(transcript, genome = gen, chromosome = chr, name = "Transcripts", cex.group = 0.5,
                           transcriptAnnotation = "transcript", fontsize.title = fontsize.title, col = "black", 
                           fill = "#85c0f9")
grtrack_sig <- HighlightTrack(trackList = grtrack, chromosome = chr, inBackground = TRUE, alpha = 1, col = "yellow",
                              fill = "yellow", start = DEU_table$start - 100, end = DEU_table$end + 100)


#==== PLOT DIFFERENTIAL HISTONE MODIFICATION ====

#---- 1. Prepare ChIP_seq signals and MAnorm results for plotting ---- 
prepare_DHM_table <- function(epi_id1, epi_id2, gene, histone_type) {
  #' Generate the datatable containing M-values and significant DHM events for each histone context
  #' 
  #' This function parse the files where M-values are annotated to exon or intron region 
  #' @param epi_id1 name of the first tissue in the comparison
  #' @param epi_id1 name of the second tissue in the comparison
  #' @param gene gene name to plot
  #' @param histone_type histone context of the comparison
  #' @return an GRanges object containing annotated M-values and p-values of detected ChIP-seq peaks
  #' 
  
  # Parse files for M-values annotation in exon and intron regions
  file_id_exon <- paste(".*exon_", paste(histone_type, "txt", sep = "."), sep = ".*")
  file_id_pi <- paste(".*pi_", paste(histone_type, "txt", sep = "."), sep = ".*")
  DHM_exon <- fread(normalizePath(all_files[grep(file_id_exon, all_files, fixed = FALSE)]),
                    header = FALSE, stringsAsFactors = FALSE, quote = FALSE)
  DHM_pi <- fread(normalizePath(all_files[grep(file_id_pi, all_files, fixed = FALSE)]),
                  header = FALSE, stringsAsFactors = FALSE, quote = FALSE)
  DHM_tables_list <- list(DHM_exon, DHM_pi)

  # Get location of peaks, M-values and p-values for plotting
  for (i in seq(length(DHM_tables_list))) {
    DHM_table <- DHM_tables_list[[i]]
    colnames(DHM_table) = c("chr", "source", "feature", "start", "end", "score", "strand", "frame", "metadata",
                            "chr2", "start2", "end2", "m_val", "p_val", "signal_1", "signal_2")
    DHM_table <- DHM_table %>%
      dplyr::filter(feature == "exonic_part" | feature == "promoter" | feature == "intron") %>%
      dplyr::mutate(
        start = dplyr::if_else(start2 == -1, start, start2), end = dplyr::if_else(end2 == -1, end, end2),
        signal_1 = as.numeric(as.character(signal_1)), signal_2 = as.numeric(as.character(signal_2)),
        signal_1 = dplyr::if_else(is.na(signal_1), 0, signal_1), signal_2 = dplyr::if_else(is.na(signal_2), 0, signal_2),
        m_val = as.numeric(as.character(m_val)), p_val = as.numeric(as.character(p_val)),
        m_val = dplyr::if_else(is.na(m_val), 0, m_val), p_val = dplyr::if_else(is.na(p_val), 1, p_val)
      ) %>%
      dplyr::select(chr, feature, strand, start, end, m_val, p_val, signal_1, signal_2)
    colnames(DHM_table)[(ncol(DHM_table) - 1) : ncol(DHM_table)] = c(epi_id1, epi_id2)

    DHM_tables_list[[i]] <- DHM_table
  }
  DHM_tables_list <- rbind(DHM_tables_list[[1]], DHM_tables_list[[2]])
  DHM_tables_list <- as(DHM_tables_list, "GRanges")
  return(DHM_tables_list)
}

get_all_DHM_tables <- function(m_sig = 1, p_sig = 0.05) {
  #' This function prepare the ChIP-seq signals from each tissue, M-values, p-values and mark the DHM for plotting
  #' 
  #' @param m_sig threshold for M-values to define DHM
  #' @param p_sig threshold for p-values to define DHM 
  #' @return a list of datatables storing ChIP-seq signal and a list of datatables storing significant DHM exons for 
  #' each histone context
  #

  all_histone_counts <- vector("list", length(histone_types))
  all_histone_sigs <- vector("list", length(histone_types))
  for (i in seq(length(histone_types))) {
     # Get table containing signals for each tissue, M-values, and p_values
    histone_type <- histone_types[i]
    DHM_table <- prepare_DHM_table(epi_id1, epi_id2, gene_name, histone_type)

    DHM_table_count <- DHM_table[, c(epi_id1, epi_id2)]
    # Define significant DHMs
    DHM_table_sig <- as.data.frame((DHM_table[DHM_table$p_val <= 0.05 & abs(DHM_table$m_val) >= m_sig]))

    all_histone_counts[[i]] <- DHM_table_count
    all_histone_sigs[[i]] <- DHM_table_sig
  }
  
  return(list(all_histone_counts, all_histone_sigs))
}

all_histone_counts <- get_all_DHM_tables()
all_histone_sigs <- all_histone_counts[[2]]
all_histone_counts <- all_histone_counts[[1]]


#---- 2. Plot ChIP-seq signals with highlighted DHM exons ----  
get_all_histone_track_sigs <- function(all_histone_counts, all_histone_sigs, expansion_range) {
  #' Prepare GViz tracks for plotting from ChIP-seq data and MAnorm results 
  #' 
  #' Given dataframes 
  #' @param all_histone_counts a list of dataframes with ChIP-seq signals and their coordinates for each histone context
  #' @param all_histone_sigs a list of dataframes with M-values and p-values from MAnorm and their coordinates for each 
  #' histone context
  #' @param expansion_range the distance in bp to expand highlighting track for better visualization
  #' @return a list of GViz tracks for plotting
  #' 
  
  all_histone_track_sigs <- vector("list", length(all_histone_counts))
  # For each histone type
  for (i in seq(length(all_histone_counts))) {
    histone_count <- all_histone_counts[[i]]
    histone_sig <- all_histone_sigs[[i]]

    # If plotting the last histone, include legend
    if (i == length(all_histone_counts)) {
      histoneTrack <- DataTrack(histone_count, name = histone_types[i], type = c("a"), legend = TRUE, 
                          groups = c(epi_name1, epi_name2), col = type.color)
    } else {
      histoneTrack <- DataTrack(histone_count, name = histone_types[i], type = c("a"), legend = FALSE,
                          groups = c(epi_name1, epi_name2), col = type.color)
    }

    # If there is at least one significant DHM, plot DataTrack and HighlightTrack
    if (dim(histone_sig)[[1]] != 0) {
      histoneTrack_sig <- HighlightTrack(trackList = histoneTrack, chromosome = chr, col = col.highlight,
                                         fill = col.highlight, start = histone_sig$start - expansion_range,
                                         end = histone_sig$end + expansion_range)
      all_histone_track_sigs[[i]] <- histoneTrack_sig
    } else {
      all_histone_track_sigs[[i]] <- histoneTrack
    }
  }
  return(all_histone_track_sigs)
}
all_histone_track_sigs <- get_all_histone_track_sigs(expansion_range = 0, all_histone_counts = all_histone_counts, 
                                                     all_histone_sigs = all_histone_sigs)

#==== PLOT ALL TRACKS FOR FIGURE 4 ====
all_tracks <- list(itrack, gtrack, grtrack_sig)
all_tracks <- append(all_tracks, all_histone_track_sigs)

tiff(filename = paste(paste(epi_id1, epi_id2, gene_name, sep = "_"), "tiff", sep = "."), width = 12, height = 15, 
     units = "in", res = 300)
plotTracks(all_tracks, extend.left = 5000, extend.right = 500, fontsize.title = fontsize.title, stackHeight = 0.75,
           background.title = background.title, col.title = col.title, margin = margin, reverseStrand = reverseStrand,
           stackHeight = 0.75, cex = 0.6, lwd = 1.5, frame = TRUE, cex.axis = 0.3, lineheight = 0.25, 
           cex.legend = 0.75, fontsize.legend = fontsize.title)
dev.off()
