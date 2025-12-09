#' Annotate regions according to a provided bed annotation file 
#'
#' @import rhdf5
#' @import igraph
#' @import leidenAlg
#' @import GenomicRanges
#' @import rtracklayer
#' @import liftOver
#' @import leidenAlg
#' @import dplyr
#' @importFrom tibble tibble
#' @importFrom tidyr pivot_longer
#' 
#' @param object TChroNetNetwork object
#' @param bed_path Path to a bed file containing regions and annotations
#' @return A TChroNetNetwork object
#' @export
annotate_regions_from_bed <- function(object, bed_path , genome = 'hg39' ) {
  if (!inherits(object, "TCrhoNetNetwork")) {
    stop("The 'object' must be a TCrhoNetNetwork instance.")
  }
  if (!file.exists(bed_path)) {
    stop("Provided file does not exist.")
  }
  if (!genome %in% c("hg19" , "hg38" , "mm10")) {
    stop("Provided correct 'genome'.")
  }
  
  if (genome %in% c("hg19" , "hg38")){
    if (genome == 'hg19'){
      regions_list <- as.data.frame(object@genomicRegions)
  
      gr_network <- GenomicRanges::GRanges(
        seqnames = regions_list$seqnames  ,
        ranges = IRanges::IRanges(start = regions_list$start, end = regions_list$end),
        node = regions_list$nodes
      )
      
      annotation_results <- data.frame(node = regions_list$nodes, annotation = NA_character_, stringsAsFactors = FALSE)
    }
    else {
      regions_list <- object@lifted_coords[!is.na(object@lifted_coords$chr_hg38),]
  
      gr_network <- GenomicRanges::GRanges(
        seqnames = regions_list$chr_hg38  ,
        ranges = IRanges::IRanges(start = regions_list$start_hg38, end = regions_list$end_hg38),
        node = regions_list$lifted_coord
      )
      
      
      annotation_results <- data.frame(node = regions_list$lifted_coord, annotation = NA_character_, stringsAsFactors = FALSE)
    }
    
  }
  else {
    regions_list <- as.data.frame(object@genomicRegions)
  
    gr_network <- GenomicRanges::GRanges(
      seqnames = regions_list$seqnames  ,
      ranges = IRanges::IRanges(start = regions_list$start, end = regions_list$end),
      node = regions_list$nodes
    )
      
    annotation_results <- data.frame(node = regions_list$nodes, annotation = NA_character_, stringsAsFactors = FALSE)
  }

  
  
  bed_gr <- rtracklayer::import(bed_path, format = "BED")
  overlaps <- findOverlaps(gr_network, bed_gr)
  
  overlaps <- overlaps[!duplicated(queryHits(overlaps)),]
  
  overlapping_nodes <- gr_network$node[queryHits(overlaps)]
  resulted_annot <- bed_gr$name[subjectHits(overlaps)]
  annotation_results$annotation[annotation_results$node %in% overlapping_nodes] <- resulted_annot

  annotation_results <- annotation_results[!is.na(annotation_results$annotation),]
  # --- 3. Add to object ---
  object@annotations <- annotation_results
  
  message("✅ Annotation complete: ", length(unique(na.omit(annotation_results$annotation))), " annotation types added.")
  
  return(object)
}