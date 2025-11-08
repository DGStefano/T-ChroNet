#' Evalaute the annotations for each investigated region
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
#' @param resolutions Set the resolution parameter for plotting trends. If not specified the function the resolution parameter stored in the object (deafault = NULL)
#' @return A dataframe containing annotations for each regions
#' @export

genomic_position_stackbar <- function(object , resolution = NULL) {
  if (!inherits(object, "TCrhoNetNetwork")) {
    stop("The 'object' must be a TCrhoNetNetwork instance.")
  }
  if (is.null(resolution)){
    resolution = object@resolution
  }
  # --- 1. Identify the cluster column name ---
  cluster_col <- paste0("clusters_", as.character(resolution))
  
  if (!(cluster_col %in% colnames(object@clusters))) {
    stop(paste0("No clustering results found for resolution ", resolution))
  }

  # Initialize empty list for results
  combined_list <- list()

  # Iterate over membership_nodes
  for (i in unique(object@clusters[[cluster_col]])) {
    community_nodes <- object@clusters[object@clusters[cluster_col] == i, 'node']

    community_nodes_lifted <- object@lifted_coords[object@lifted_coords$node %in% community_nodes, "lifted_coord"]
    # Filter annotation_df for nodes in this community
    community_annot <- object@annotations[object@annotations$node %in% community_nodes_lifted, , drop = FALSE]

    # Remove NA annotations
    community_annot <- community_annot[!is.na(community_annot$annotation), , drop = FALSE]

    if (nrow(community_annot) > 0) {
      community_annot$community_number <- paste0("community_", i)
      combined_list[[i]] <- community_annot
      }
    }

  # Combine all membership_nodes into one data frame
  final_annotation_df <- do.call(rbind, combined_list)

  # Return result
  message("✅ Combined annotations from ", length(combined_list), " membership_nodes.")
  return(final_annotation_df)
} 
