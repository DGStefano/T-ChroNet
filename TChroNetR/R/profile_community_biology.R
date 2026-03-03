#' Profile Community Biology
#' 
#' @importFrom GenomicRanges GRanges promoters distanceToNearest seqnames start end width
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors metadata mcols
#' @importFrom BiocGenerics score
#' @importFrom dplyr select everything
#' 
#' @param object TChroNetNetwork object
#' @param genome A BSgenome object (e.g., BSgenome.Hsapiens.UCSC.hg38).
#' @param phast A GScores or GenomicScores object (e.g., phastCons100way.UCSC.hg38).
#' @param txdb A TxDb object for TSS calculation (e.g., TxDb.Hsapiens.UCSC.hg38.knownGene).
#' @param resolution The resolution value used for community detection
#' @return A dataframe containing statistics for each genomic region
#' @export
profile_community_biology <- function(object, genome, phast, txdb, resolution) {
  if (!inherits(object, "TCrhoNetNetwork")) {
    stop("Input must be a 'TCrhoNetNetwork' object.")
  }

  if (missing(resolution)) {
    if (length(object@resolution) > 0) {
      resolution <- object@resolution
    } else {
      stop("Resolution parameter is missing and not found in the object.")
    }
  }

  is_hg38 <- grepl("hg38", metadata(genome)$genome, ignore.case = TRUE)
  
  if (is_hg38) {
    if (nrow(object@lifted_coords) == 0) {
      stop("Genome is hg38 but 'lifted_coords' slot is empty.")
    }
    
    valid_indices <- which(!is.na(object@lifted_coords$lifted_coord))
    
    if (length(valid_indices) == 0) {
      stop("All lifted_coords are NA.")
    }
    
    lifted_subset <- object@lifted_coords[valid_indices, ]
    
    gr <- GRanges(
      seqnames = lifted_subset$chr_hg38,
      ranges = IRanges(
        start = as.numeric(lifted_subset$start_hg38),
        end = as.numeric(lifted_subset$end_hg38)
      )
    )
    
    df_total <- data.frame(node = lifted_subset$node)
    
  } else {
    gr <- object@genomicRegions
    df_total <- data.frame(node = names(gr))
    if(is.null(df_total$node)) df_total$node <- seq_along(gr)
  }

  all_tss <- GenomicRanges::promoters(txdb, upstream = 0, downstream = 1)
  
  df_total$chr <- as.character(seqnames(gr))
  df_total$start <- start(gr)
  df_total$end <- end(gr)
  df_total$region_length <- GenomicRanges::width(gr)
  
  seqs <- BSgenome::getSeq(genome, gr)
  df_total$gc_content <- as.numeric(Biostrings::letterFrequency(seqs, letters = "GC", as.prob = TRUE))
  
  dist_to_tss <- GenomicRanges::distanceToNearest(gr, all_tss)
  df_total$tss_distance <- mcols(dist_to_tss)$distance
  
  df_total$score_phylo <- BiocGenerics::score(phast, gr)
  
  community_col <- paste0("clusters_", as.character(resolution))
  
  if (!community_col %in% names(object@clusters)) {
    stop(paste("Resolution", resolution, "not found in object@clusters."))
  }
  
  full_clusters <- object@clusters[[community_col]]
  
  if (!is.null(names(full_clusters))) {
    df_total$community_id <- full_clusters[as.character(df_total$node)]
  } else {
    df_total$community_id <- full_clusters[valid_indices]
  }

  return(df_total %>% select(node, chr, start, end, community_id, everything()))
}