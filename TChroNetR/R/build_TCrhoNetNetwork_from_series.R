#' Convert a TCrhoNetSeries into a TCrhoNetNetwork 
#' @import igraph
#' @import GenomicRanges
#' @import rtracklayer
#' @import rhdf5
#' @import dplyr
#' @importFrom data.table fread
#' @importFrom S4Vectors IRanges
#' @importFrom tibble as_tibble
#' @importFrom tidyr pivot_longer
#' 
#' @param series_object TCrhoNetSeries object
#' @param threshold Choosen correlation threshold. If not provided it uses the best threshold parameter stored in the TCrhoNetSeries object
#' @param matrix_path Path of the normalized log2 matrix used for correlation evaluation
#' @param verbose (default=TRUE)
#' @return A TCrhoNetNetwork object
#' @export
#' Optimized conversion of TCrhoNetSeries into a TCrhoNetNetwork 
build_TCrhoNetNetwork_from_series <- function(series_object,
                                              threshold,
                                              matrix_path,
                                              verbose = TRUE) {
  
  # --- Input Validation & Thresholding ---
  if (!inherits(series_object, "TCrhoNetSeries")) {
    stop("❌ Input must be a 'TCrhoNetSeries' object.")
  }
  
  if (missing(threshold)) {
    if ('best_th' %in% slotNames(series_object)) {
      threshold <- series_object@best_th
    } else {
      stop("❌ You must specify a 'threshold' or run 'find_best_th' first.")
    }
  }

  # --- Graph Filtering ---
  if (verbose) message("⚡ Filtering graph at threshold ", threshold, "...")
  
  # --- Subgraph to remove edges below the threshold ---
  edge_mask <- E(series_object@graph)$weight >= threshold
  G <- igraph::subgraph_from_edges(series_object@graph, which(edge_mask), delete.vertices = FALSE)

  # --- Reading matrix ---
  if ('matrix' %in% slotNames(series_object)) {
    matrix_G <- series_object@matrix
  } else {
    if (verbose) message("⚡ Loading matrix using ...")
    raw_data <- data.table::fread(matrix_path, data.table = FALSE, sep = "\t")
    
    # Set rownames and remove the first column
    rownames(matrix_G) <- raw_data[[1]]
    matrix_G <- as.matrix(raw_data[, -1]) 
    rm(raw_data) # Free memory
  }

  # --- Creation of genomic choordinates---
  if (verbose) message("⚡ Parsing genomic coordinates...")
  
  node_names <- rownames(matrix_G)
  
  coords <- do.call(rbind, strsplit(node_names, "-"))
  
  peaks_ranges <- GenomicRanges::GRanges(
    seqnames = coords[, 1],
    ranges = IRanges::IRanges(
      start = as.numeric(coords[, 2]),
      end = as.numeric(coords[, 3])
    ),
    nodes = node_names
  )

  # --- Object Assembly ---
  th_str <- as.character(threshold)
  
  network_obj <- new(
    "TCrhoNetNetwork",
    graph = G,
    clusters = if(th_str %in% names(series_object@communities)) series_object@communities[[th_str]] else list(),
    matrix = matrix_G,
    genomicRegions = peaks_ranges,
    modularity = if(th_str %in% names(series_object@modularity)) as.data.frame(series_object@modularity[[th_str]]) else data.frame(),
    metadata = list(
      n_nodes = igraph::vcount(G), 
      n_edges = igraph::ecount(G),
      threshold = threshold
    )
  )

  if (verbose) {
    message(sprintf(
      "✅ Success: %d nodes, %d edges preserved.",
      igraph::vcount(G), igraph::ecount(G)
    ))
  }
  
  return(network_obj)
}