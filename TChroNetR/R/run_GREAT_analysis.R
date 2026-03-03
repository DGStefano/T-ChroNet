#' Find taret genes for each community accordig GREAT analysis
#' 
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom GenomeInfoDb seqnames
#' 
#' @param object TChroNetNetwork object
#' @param resolutions Set the resolution parameter for plotting trends. If not specified the function the resolution parameter stored in the object (deafault = NULL)
#' @param genome Genome to investigate (default = 'hg19')
#' @return A TChroNetNetwork object
#' @export

run_GREAT_analysis <- function(object,
                               resolution = NULL,
                               genome = "hg19"
                              ) {
  
  if (!inherits(object, "TCrhoNetNetwork")) {
    stop("The 'object' must be a TCrhoNetNetwork instance.")
  }
  
  if (is.null(object@clusters) || !"data.frame" %in% class(object@clusters)) {
    stop("The object must contain a valid @clusters data frame.")
  }
  
  if (is.null(object@genomicRegions) && is.null(object@lifted_coords)) {
    stop("The object must have genomic coordinates in @genomicRegions or @lifted_coords.")
  }
  
  suppressPackageStartupMessages({
    library(rGREAT)
    library(GenomicRanges)
    library(dplyr)
  })

  if(is.null(resolution)) {
    resolution <- object@resolution
  }
  
  cluster_col <- paste0('clusters_' , as.character(resolution))
  
  message("🔍 Using cluster column: ", cluster_col)
  
  # --- 2. Prepare coordinates ---
  if (!is.null(object@lifted_coords) & nrow(object@lifted_coords) > 0 ) {
    coords_df <- object@lifted_coords
  } else {
    coords_df <- data.frame(
      chr = as.character(GenomeInfoDb::seqnames(object@genomicRegions)),
      start = start(object@genomicRegions),
      end = end(object@genomicRegions),
      node = object@genomicRegions$nodes
    )
  }
  
  clusters_df <- object@clusters
  clusters_df <- clusters_df[, c("node", cluster_col)]
  colnames(clusters_df) <- c("node", "cluster")
  
  # --- 3. Merge coordinates and cluster assignments ---
  merged_df <- merge(coords_df, clusters_df, by = "node", all.x = TRUE)
  merged_df <- merged_df[!is.na(merged_df$cluster), ]
  
  unique_membership_nodes <- sort(unique(merged_df$cluster))
  message("📊 Found ", length(unique_membership_nodes), " membership_nodes to analyze.")
  
  # --- 4. Initialize storage ---
  final_target_genes <- list()
  
  # --- 5. Loop over membership_nodes ---
  for (comm in unique_membership_nodes) {
    message("→ Running GREAT for community ", comm, " ...")
    
    comm_regions <- merged_df[merged_df$cluster == comm, ]
    if (nrow(comm_regions) == 0) next
    
    gr <- makeGRangesFromDataFrame(
      comm_regions,
      seqnames.field = colnames(comm_regions)[2],
      start.field = colnames(comm_regions)[3],
      end.field = colnames(comm_regions)[4],
      keep.extra.columns = TRUE
    )
    
    GREAT_results <- interrogate_GREAT_regions(
      gr = gr,
      comm_name = paste0("community_", comm),
      species = genome
    )
    
   
    final_target_genes[[comm]] <- GREAT_results
  }
  
  # --- 6. Attach to object ---
  object@GREAT_targets <- final_target_genes
  
  message("✅ GREAT analysis complete: ", length(final_target_genes), " membership_nodes processed.")
  return(object)
}