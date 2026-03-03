#' Create the TChroNetSeries object
#' @importFrom igraph make_empty_graph add_vertices add_edges vcount ecount E V degree components edge_density transitivity
#' @importFrom leidenbase leiden_find_partition
#' @importFrom dplyr select collect mutate filter arrange left_join %>%
#' @importFrom arrow open_dataset
#' @importFrom methods new
#' @importFrom stats setNames
#' 
#' @param output_dir Directory containing TChroNetPy results
#' @param matrix_path Path to the log2 normalized count matrix used in tchronetPy
#' @param method Method for community investigation (c("Leiden", "Louvain"))
#' @param resolutions_list List of resolutions to be investigated 
#' @param run_cd Boolean value for running community detection evaluation. default = FALSE
#' @param transitivity Boolean value for evaluating transitivity value for each threshold. default = FALSE
#' @param min_size Size of the smallest community to consider
#' @param merged_id The community id given to all the nodes belonging to communities smaller than min_size
#' @param seed Set a seed for computational reproducibility
#' @return A TChroNetSeries object
#' @export
build_TChroNetSeries_object <- function(output_dir, matrix_path,
                                        method = c("Leiden", "Louvain"),
                                        resolutions_list = seq(0.1, 0.9, 0.1),
                                        run_cd = FALSE,
                                        transitivity = FALSE,
                                        min_size = 100,
                                        merged_id = 1000,
                                        seed = 123) {
  # --- Required Libraries ---
  library(igraph)
  library(dplyr)
  library(arrow)
  library(leidenbase)

  # --- Internal Helper: Fix Leiden Membership ---
  .fix_leiden_membership <- function(G, membership, min_size = 100, merged_id = 1000) {
    if (length(membership) < igraph::vcount(G)) {
      full <- rep(NA_integer_, igraph::vcount(G))
      full[seq_along(membership)] <- membership
      membership <- full
    }
    # Isolated nodes get unique singleton IDs
    iso <- which(igraph::degree(G) == 0)
    if (length(iso) > 0) {
      max_id <- max(membership, na.rm = TRUE)
      if (is.infinite(max_id)) max_id <- 0
      membership[iso] <- seq_len(length(iso)) + max_id
    }
    # Handle small communities
    comm_sizes <- table(membership)
    small_comms <- as.integer(names(comm_sizes[comm_sizes < min_size]))
    if (length(small_comms) > 0) {
      membership[membership %in% small_comms] <- merged_id
    }
    return(membership)
  }

  method <- match.arg(method)
  if (!run_cd) method <- ""

  # --- Scan and Sort Bin Directories ---
  bin_dirs <- list.dirs(output_dir, full.names = TRUE, recursive = FALSE)
  bin_dirs <- bin_dirs[grepl("bin_", bin_dirs)]
  if (length(bin_dirs) == 0) stop("No 'bin_' folders found.")

  extract_thr <- function(p) as.numeric(sub(".*bin_(\\d+\\.\\d+)", "\\1", basename(p)))
  bin_dirs <- bin_dirs[order(-extract_thr(bin_dirs))]
  thresholds <- extract_thr(bin_dirs)

  # --- Load Matrix with EXACT Names ---
  message("📂 Loading matrix (check.names = FALSE)...")
  matrix_G <- read.delim(matrix_path, sep = "\t", row.names = 1, check.names = FALSE)
  node_names <- as.character(rownames(matrix_G))

  # --- Initialize Object ---
  obj <- new("TCrhoNetSeries")
  obj@method <- method
  obj@matrix <- matrix_G
  obj@communities <- list()
  obj@modularity <- list()
  obj@metrics <- data.frame()

  # Initialize graph with all nodes
  G <- igraph::make_empty_graph(directed = FALSE)
  G <- igraph::add_vertices(G, length(node_names), name = node_names)

  # --- Iterative Build ---
  for (i in seq_along(bin_dirs)) {
    thr_path <- bin_dirs[i]
    thr <- thresholds[i]
    message("\n🔹 Threshold Bin: ", thr)

    # Robust Arrow Load
    e_df <- tryCatch({
      ds <- arrow::open_dataset(thr_path, format = "parquet")
      df <- ds %>% 
        dplyr::select(row_name_s, row_name_t, corr) %>% 
        dplyr::collect() %>%
        dplyr::mutate(
          row_name_s = trimws(as.character(row_name_s)),
          row_name_t = trimws(as.character(row_name_t))
        )
      
      # IDENTITY FIX 1: Filter to ensure nodes exist in matrix
      df <- df %>% filter(row_name_s %in% node_names & row_name_t %in% node_names)
      
      # IDENTITY FIX 2: Deterministic sort for identical Leiden results
      # This ensures the graph adjacency list is built in the exact same order
      df <- df %>% dplyr::arrange(row_name_s, row_name_t)
      
      df
    }, error = function(e) {
      message("❌ Read Error: ", e$message)
      return(NULL)
    })

    if (is.null(e_df) || nrow(e_df) == 0) {
      message("   (No valid edges in bin)")
      next
    }

    # Add edges to the cumulative graph
    edge_vector <- as.vector(t(as.matrix(e_df[, c("row_name_s", "row_name_t")])))
    G <- igraph::add_edges(G, edge_vector, attr = list(weight = e_df$corr))

    if (igraph::ecount(G) == 0) next

    # --- Metrics ---
    comp <- igraph::components(G, mode = "weak")
    lcc_size <- max(comp$csize)
    n_nodes <- igraph::vcount(G)
    relative_lcc <- lcc_size / n_nodes

    trans_val <- 0
    if (transitivity) {
      trans_val <- suppressWarnings(igraph::transitivity(G, type = "global"))
    }

    # --- Community Detection ---
    th_cluster_df <- data.frame(node = igraph::V(G)$name, stringsAsFactors = FALSE)
    modularity_at_th <- data.frame()

    if (method == "Leiden") {
      for (res in resolutions_list) {
        message(sprintf("   🧩 Leiden [res: %.2f]...", res))
        
        partition <- leidenbase::leiden_find_partition(
          igraph = G, 
          edge_weights = igraph::E(G)$weight,
          partition_type = "RBConfigurationVertexPartition",
          resolution_parameter = res, 
          seed = seed
        )

        membership_nodes <- .fix_leiden_membership(G, partition$membership, min_size, merged_id)
        
        cluster_col <- paste0("clusters_", res)
        tmp_df <- data.frame(node = igraph::V(G)$name, clusters = membership_nodes)
        colnames(tmp_df)[2] <- cluster_col
        
        th_cluster_df <- dplyr::left_join(th_cluster_df, tmp_df, by = "node")
        
        modularity_at_th <- rbind(modularity_at_th, data.frame(
          resolution = res, 
          modularity = partition$modularity, 
          quality = partition$quality, 
          threshold = thr,
          stringsAsFactors = FALSE
        ))
      }
    }
    if(run_cd == FALSE) {
      membership_nodes = NA
    }
    # --- Save Results ---
    obj@communities[[as.character(thr)]] <- th_cluster_df
    obj@modularity[[as.character(thr)]] <- modularity_at_th
    obj@metrics <- rbind(obj@metrics, data.frame(
      threshold = thr, 
      n_nodes = igraph::vcount(G), 
      n_edges = igraph::ecount(G),
      density = igraph::edge_density(G),
      transitivity = trans_val,
      lcc = lcc_size, 
      relative_lcc = relative_lcc,
      n_communities = ifelse(all(is.na(membership_nodes)), 0, length(unique(membership_nodes))),
      stringsAsFactors = FALSE
    ))
    
    message(sprintf("   ✅ Added %d edges. Total edges: %d", nrow(e_df), igraph::ecount(G)))
  }
  
  obj@graph <- G
  message("\n🏁 TChroNet series construction complete.")
  return(obj)
}