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
compute_membership_nodes <- function(object, resolutions = 1, method = c("Leiden", "Louvain") , niter = 2 , min_size = 100, merged_id = 1000, seed = 1234) {

  .add_weighted_self_loops_iso <- function(G, loop_weight = 1) {
    iso <- which(degree(G) == 0)
    if (length(iso) == 0) return(G)
    add_edges(
      G,
      as.vector(rbind(iso, iso)),
      attr = list(weight = rep(loop_weight, length(iso)))
    )
  }

  .fix_leiden_membership <- function(G, membership, min_size = 100, merged_id = 1000) {

    # Ensure full length
    if (length(membership) < vcount(G)) {
      full <- rep(NA_integer_, vcount(G))
      full[seq_along(membership)] <- membership
      membership <- full
    }

    # Handle degree-0 nodes → unique singleton communities
    iso <- which(degree(G) == 0)
    if (length(iso) > 0) {
      max_id <- max(membership, na.rm = TRUE)
      membership[iso] <- seq_len(length(iso)) + max_id
    }

    # Community sizes
    comm_sizes <- table(membership)

    # Communities smaller than min_size
    small_comms <- as.integer(names(comm_sizes[comm_sizes < min_size]))

    if (length(small_comms) > 0) {
      membership[membership %in% small_comms] <- merged_id
    }

    return(membership)
  }
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
  
    G <- object@graph
  
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
          

        partition_current_graph <- leidenbase::leiden_find_partition(
          igraph = G,
          edge_weights = E(G)$weight,
          partition_type = "RBConfigurationVertexPartition", # Options: ModularityVertexPartition, etc.
          resolution_parameter = res,
          seed = seed
        )
                
        membership_nodes_print <- .fix_leiden_membership(
          G,
          membership_nodes,
          min_size = min_size,
          merged_id = merged_id
        )

        cluster_col <- paste0("clusters_", res)
        tmp_df <- data.frame(node = V(G)$name, stringsAsFactors = FALSE)
        tmp_df[[cluster_col]] <- membership_nodes_print
          
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
          mod_val <- partition_current_graph$modularity
          qual_val <- partition_current_graph$quality
          object@modularity <- rbind(object@modularity, data.frame(
            resolution = res,
            modularity = mod_val,
            quality = qual_val,
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
        resolution = res,
        modularity = mod_val,
        method = "Louvain",
        stringsAsFactors = FALSE
      ))
    }
    return(object)
}