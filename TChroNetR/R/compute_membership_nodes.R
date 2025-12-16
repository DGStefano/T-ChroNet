#' Calcualte communities in network
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
#' @param resolutions Resolution parameter it can be a unique resolution or a list of resolutions
#' @param method Choose a 'Leiden' or 'Louvain' method for community detection
#' @param niter Number of Leiden iterations that the algorithm should be run for (default=2)
#' @param seed Set a seed for computational reproducibility
#' @return An object of class TChroNetNetwork
#' @export
compute_membership_nodes <- function(object, resolutions = 1, method = c("Leiden", "Louvain") , niter = 2 , seed = 1234) {

    method <- match.arg(method)
    
    if (!inherits(object@graph, "igraph")) {
      stop("The 'graph' slot must contain a valid igraph object.")
    }
    
    # Ensure vertex names exist
    if (is.null(V(object@graph)$name)) {
      V(object@graph)$name <- as.character(seq_len(vcount(object@graph)))
    }
    
    # Ensure edge weights are numeric
    if (!is.null(E(object@graph)$corr)) {
      E(object@graph)$corr <- as.numeric(E(object@graph)$corr)
    }
    
    # Initialize modularity slot if missing
    if (is.null(slot(object, "modularity")) || 
        !is.data.frame(object@modularity)) {
      object@modularity <- data.frame()
    }
    
    # Initialize clusters slot if missing
    if (is.null(slot(object, "clusters")) || 
        !is.data.frame(object@clusters)) {
      object@clusters <- data.frame()
    }
    set.seed(seed)
    # --- LEIDEN METHOD ---
    if (method == "Leiden") {
      for (res in resolutions) {
        if (paste0("clusters_", as.character(res)) %in% colnames(object@clusters)) {
          message(sprintf("Resolution %.2f already calculated ", res))
          next
        }
        else {
          message(sprintf("Running Leiden clustering at resolution %.2f...", res))
          
          set.seed(seed)
 
          nodes_membership <- leidenAlg::find_partition_with_rep(object@graph, resolution = res, edge_weights = E(object@graph)$weight , nrep = 5)

          nodes_membership <- nodes_membership+1

          if (length(nodes_membership$membership) != vcount(object@graph)) {
            stop(sprintf(
              "Leiden returned %d assignments but the graph has %d vertices.",
              length(nodes_membership$membership), vcount(object@graph)
            ))
          }
          # Build data frame for this resolution
          cluster_col <- paste0("clusters_", as.character(res))
          tmp_df <- data.frame(
            node = V(object@graph)$name,
            stringsAsFactors = FALSE
          )
          tmp_df[[cluster_col]] <- as.integer(nodes_membership$membership)
          
          # Handle first run (clusters empty)
          if (nrow(object@clusters) == 0) {
            object@clusters <- tmp_df
          } else if (!"node" %in% colnames(object@clusters)) {
            # if clusters exist but no node column (edge case)
            object@clusters$node <- V(object@graph)$name
            object@clusters <- merge(object@clusters, tmp_df, by = "node", all = TRUE)
          } else {
            # Safe merge by node
            object@clusters <- merge(object@clusters, tmp_df, by = "node", all = TRUE)
          }
          
          # Compute modularity
          mod_val <- igraph::modularity(object@graph, membership = nodes_membership$membership, weights = E(object@graph)$corr , directed = FALSE)
          object@modularity <- rbind(object@modularity, data.frame(
            resolution = res,
            modularity = mod_val,
            method = "Leiden",
            stringsAsFactors = FALSE
          ))
          object@resolution <- res
        }
      }
    }
    
    # --- LOUVAIN METHOD ---
    if (method == "Louvain") {
      message("Running Louvain clustering...")
      
      membership_nodes <- cluster_louvain(object@graph, weights = E(object@graph)$corr)
      membership <- membership(membership_nodes)
      
      tmp_df <- data.frame(
        node = names(membership),
        clusters_Louvain = as.integer(membership),
        stringsAsFactors = FALSE
      )
      
      if (nrow(object@clusters) == 0) {
        object@clusters <- tmp_df
      } else {
        object@clusters <- merge(object@clusters, tmp_df, by = "node", all = TRUE)
      }
      
      mod_val <- modularity(membership_nodes)
      object@modularity <- rbind(object@modularity, data.frame(
        resolution = 1,
        modularity = mod_val,
        method = "Louvain",
        stringsAsFactors = FALSE
      ))
    }
    return(object)
}