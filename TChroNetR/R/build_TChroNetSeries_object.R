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
#' @param edge_files List of parquet files to build the network
#' @param matrix_path Path to the log2 normalized count matrix used in tchronetPy
#' @param method Method for community investigation (c("Leiden", "Louvain"))
#' @param resolutions_list List of resolutions to be investigated 
#' @param run_cd Boolean value for running community detection evaluation. default = FALSE
#' @param transitivity Boolean value for evaluating transitivity value for each threshold. default = FALSE
#' @param seed Set a seed for computational reproducibility
#' @return A TChroNetSeries object
#' @export
#' 

build_TChroNetSeries_object <-  function(edge_files, matrix_path,
                                      method = c("Leiden", "Louvain"),
                                      resolutions_list = seq(0.1, 0.9, 0.1),
                                      run_cd = FALSE,
                                      transitivity = FALSE,
                                      min_size = 100,
                                      merged_id = 1000,
                                      seed = 123) {
  library(igraph)
  library(leidenAlg)
  library(dplyr)
  library(arrow)

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
  if (!run_cd) method <- ""


  # --- Validate inputs ---
  if (length(edge_files) == 0) stop("No edge files provided.")
  if (!all(file.exists(edge_files))) stop("Some files do not exist.")

  message("📂 Processing ", length(edge_files), " edge files...")

  # --- Sort files by descending threshold ---
  extract_thr <- function(fname) as.numeric(sub(".*_(\\d+\\.\\d+)_.*", "\\1", basename(fname)))
  edge_files <- edge_files[order(-extract_thr(edge_files))]
  thresholds <- gsub(".*edges_(\\d+\\.\\d+)_.*\\.parquet", "\\1", basename(edge_files))
  thresholds <- as.numeric(thresholds)
  trans_val = 0
  # --- Load original matrix log2 ---
  matrix_G <- read.delim(matrix_path , sep ="\t" ,row.names=1)

  # --- Initialize object ---
  obj <- new("TCrhoNetSeries")
  obj@method <- method
  obj@matrix <- matrix_G
  obj@communities <- list()
  obj@modularity <- list()
  obj@metrics <- data.frame()


  G <- igraph::make_empty_graph(directed = FALSE)
  G <- add_vertices(G, length(rownames(matrix_G)),
                  name = rownames(matrix_G) )

  for (i in seq_along(edge_files)) {
    # Step 1
    thr <- thresholds[i]
    message("🔹 Processing threshold: ", thr)
    # Step 1
    # --- Read full Parquet file ---
    e_df <- tryCatch({
      df <- arrow::read_parquet(edge_files[i])
      # convert Arrow binary to character if needed
      if (inherits(df$row_name_s, "arrow_binary")) df$row_name_s <- as.character(df$row_name_s)
      if (inherits(df$row_name_t, "arrow_binary")) df$row_name_t <- as.character(df$row_name_t)
      df <- df %>% select(row_name_s, row_name_t, corr)
      df
    }, error = function(e) {
      message("⚠️ Could not read ", edge_files[i], ": ", e$message)
      return(NULL)
    })
    if (is.null(e_df)) next
    # --- Add current edges to cumulative graph ---
    # cum_edges <- unique(dplyr::bind_rows(cum_edges, e_df))
    # --- Build graph ---
    G <- igraph::add_edges( G , t(e_df[,c(1,2)]) , attr = list(weight = e_df$corr))
    # return(G)
    # --- Compute transitivity ---
    if (transitivity == T ) {
      trans_val <- suppressWarnings(igraph::transitivity(G, type = "global"))
    }
    # --- Compute Percolation ---
    comp <- igraph::components(G, mode = "weak")
    lcc_size = max(comp$csize)
    n_nodes <- igraph::vcount(G)
    relative_lcc <- lcc_size / n_nodes
    # --- Initialize storage ---
    th_cluster_df <- data.frame()
    modularity_at_th <- data.frame()

    if (method == "Leiden") {
      for (res in resolutions_list) {
        message(sprintf("🧩 Leiden clustering at resolution %.2f...", res))

        
        G_loop <- .add_weighted_self_loops_iso(G, loop_weight = 1)

        set.seed(seed)
        membership_nodes <- leidenAlg::find_partition_with_rep(
          G_loop,
          resolution = res,
          edge_weights = E(G_loop)$weight,
          nrep = 5
        )

        rm(G_loop)
        invisible(gc())
        
        membership_nodes <- membership_nodes+1
        
        membership_nodes_print <- .fix_leiden_membership(
          G,
          membership_nodes,
          min_size = min_size,
          merged_id = merged_id
        )

        cluster_col <- paste0("clusters_", res)
        tmp_df <- data.frame(node = V(G)$name, stringsAsFactors = FALSE)
        tmp_df[[cluster_col]] <- membership_nodes_print

        if (nrow(th_cluster_df) == 0) th_cluster_df <- tmp_df
        else th_cluster_df <- merge(th_cluster_df, tmp_df, by = "node", all = TRUE)

        mod_val <- suppressWarnings(
          igraph::modularity(G, membership = membership_nodes, weights = E(G)$weight, resolution = res)
        )

        modularity_at_th <- rbind(
          modularity_at_th,
          data.frame(resolution = res, modularity = mod_val, threshold = thr, stringsAsFactors = FALSE)
        )
      }
    }
    else {
      membership_nodes_print = NA
    }

    # Store results
    obj@communities[[as.character(thr)]] <- th_cluster_df
    obj@modularity[[as.character(thr)]] <- modularity_at_th
    obj@metrics <- rbind(
      obj@metrics,
      data.frame(
        threshold = thr,
        n_nodes = vcount(G),
        n_edges = ecount(G),
        density = igraph::edge_density(G),
        transitivity = trans_val,
        lcc = lcc_size,
        relative_lcc = relative_lcc,
        n_communities = ifelse(all(is.na(membership_nodes_print)), 0, length(unique(membership_nodes_print))),
        stringsAsFactors = FALSE
      )
    )
    obj@graph <- G

    message(sprintf("✅ %s: %d nodes, %d edges, lcc=%.1f",
                    thr, vcount(G), ecount(G), lcc_size))
  }

  message("🏁 Network series completed with ", nrow(obj@metrics), " thresholds.")
  obj@edge_files <- edge_files

  return(obj)
}