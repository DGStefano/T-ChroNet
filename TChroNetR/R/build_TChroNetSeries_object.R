#' Create the TChroNetSeries object
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
#' @param edge_files List of hd5 files to build the network
#' @param method Method for community investigation (c("Leiden", "Louvain"))
#' @param resolutions_list List of resolutions to be investigated 
#' @param run_cd Boolean value for running community detection evaluation. default = FALSE
#' @param seed Set a seed for computational reproducibility
#' @return A TChroNetSeries object
#' @export
build_TChroNetSeries_object <- function(edge_files,
                                 method = c("Leiden", "Louvain"),
                                 resolutions_list = seq(0.1, 0.9, 0.1),
                                 run_cd = FALSE,
                                 seed = 123) {
  method <- match.arg(method)
  if(!run_cd) {
    method <- ""
  }

  set.seed(seed)

  # --- Validate inputs ---
  if (length(edge_files) == 0) stop("No edge files provided.")
  if (!all(file.exists(edge_files))) stop("Some provided files do not exist.")

  message("📂 Loading ", length(edge_files), " edge files...")

  # --- Sort files by descending threshold ---
  extract_thr <- function(fname) {
    as.numeric(sub(".*_(\\d+\\.\\d+)_.*", "\\1", basename(fname)))
  }
  edge_files <- edge_files[order(-extract_thr(edge_files))]
  thresholds <- gsub(".*edges_(\\d+\\.\\d+)_.*\\.h5", "\\1", basename(edge_files))

  # --- Create new object ---
  obj <- new("TCrhoNetSeries")
  obj@method <- method
  obj@communities <- list()
  obj@modularity <- list()
  obj@metrics <- data.frame()

  # --- Initialize ---
  cum_edges <- data.frame()
  graphs_list <- list()

  for (i in seq_along(edge_files)) {
    thr <- thresholds[i]
    message("🔹 Processing threshold: ", thr)

    e <- tryCatch({
      rhdf5::h5read(edge_files[i], "/")
    }, error = function(e) {
      message("⚠️ Could not read ", edge_files[i])
      return(NULL)
    })
    if (is.null(e)) next

    # Convert to data frame with standardized columns
    e_df <- as.data.frame(e)
    e_df <- e_df[, c(2:ncol(e_df))]
    colnames(e_df) <- c("origin", "target", "corr")

    # --- Add current edges to cumulative graph ---
    cum_edges <- unique(rbind(cum_edges, e_df))

    # --- Build graph ---
    G <- graph_from_data_frame(cum_edges, directed = FALSE)
    E(G)$weight <- as.numeric(cum_edges$corr)

    # --- Compute transitivity ---
    trans_val <- suppressWarnings(transitivity(G, type = "global"))

    # --- Initialize storage ---
    th_cluster_df <- data.frame()
    modularity_at_th <- data.frame()
    components <- data.frame()
    membership_vec <- rep(NA, vcount(G))

    if (method == "Leiden") {
      for (res in resolutions_list) {
        message(sprintf("Running Leiden clustering at resolution %.2f...", res))
        membership_nodes <- suppressWarnings(
          leidenAlg::find_partition(
            G,
            resolution = res,
            edge_weights = E(G)$weight
          )
        )

        # Convert from 0-based to 1-based indexing
        membership_nodes <- as.integer(membership_nodes) + 1

        cluster_col <- paste0("clusters_", as.character(res))
        tmp_df <- data.frame(
          node = V(G)$name,
          stringsAsFactors = FALSE
        )
        tmp_df[[cluster_col]] <- membership_nodes

        # Merge cluster assignments
        if (nrow(th_cluster_df) == 0) {
          th_cluster_df <- tmp_df
        } else {
          th_cluster_df <- merge(th_cluster_df, tmp_df, by = "node", all = TRUE)
        }

        # Compute modularity
        mod_val <- suppressWarnings(
          igraph::modularity(G, membership = membership_nodes, weights = E(G)$weight , resolution = 1)
        )

        modularity_at_th <- rbind(
          modularity_at_th,
          data.frame(
            resolution = res,
            modularity = mod_val,
            threshold = thr,
            stringsAsFactors = FALSE
          )
        )
      }
      components <- rbind(
        components,
        data.frame(
          th = thr,
          components = igraph::count_components(G , mode = "weak")
        )
      )
    }

    # Store into object
    obj@communities[[thr]] <- th_cluster_df
    obj@modularity[[thr]] <- modularity_at_th
    obj@components

    # Append metrics
    obj@metrics <- rbind(
      obj@metrics,
      data.frame(
        threshold = thr,
        n_nodes = vcount(G),
        n_edges = ecount(G),
        density = igraph::edge_density(G),
        transitivity = trans_val,
        n_communities = length(unique(membership_vec)),
        stringsAsFactors = FALSE
      )
    )

    message(sprintf("✅ %s: %d nodes, %d edges, transitivity=%.3f",
                    thr, vcount(G), ecount(G), trans_val))
  }
  message("🏁 Network series completed with ", nrow(obj@metrics), " thresholds.")
  obj@edge_files <- edge_files
  return(obj)
}