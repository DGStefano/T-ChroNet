#' Convert a TCrhoNetSeries into a TCrhoNetNetwork 
#' @import rhdf5
#' @import igraph
#' @import leidenAlg
#' @import GenomicRanges
#' @import rtracklayer
#' @import liftOver
#' @import leidenAlg
#' @import dplyr
#' @importFrom tibble
#' @importFrom tidyr pivot_longer
#' 
#' @param object TCrhoNetSeries object
#' @param threshold Choosen correlation threshold. If not provided it uses the best threshold parameter stored in the TCrhoNetSeries object
#' @param matrix_path Path of the normalized log2 matrix used for correlation evaluation
#' @param verbose (default=TRUE)
#' @return A TCrhoNetNetwork object
#' @export
build_TCrhoNetNetwork_from_series <- function(series_object,
                                              threshold ,
                                              matrix_path,
                                              verbose = TRUE) {
  # --- Load required packages ---
#  require(igraph)
#  require(rhdf5)
#  require(dplyr)
  
  # --- Input checks ---
  if (!inherits(series_object, "TCrhoNetSeries")) {
    stop("❌ Input must be a 'TCrhoNetSeries' object.")
  }
  
  if (missing(threshold) & !('best_th' %in% slotNames(series_object))  ) {
    stop("❌ You must specify a 'threshold' (e.g., 0.7). or run 'find_best_th'")
  }
  if(missing(threshold) & ('best_th' %in% slotNames(series_object) ) ){
    threshold = series_object@best_th
  }
  if(missing(matrix_path) & !('matrix' %in% slotNames(series_object)) ){
    stop("❌ You must specify the matrix path")
  }
  # --- Create TCrhoNetNetwork object ---
  if (verbose) message("Creating TCrhoNetNetwork object at threshold ", threshold, "...")
  
  G <- delete_edges(series_object@graph, E(series_object@graph)[weight < threshold])
  # G <- delete_vertices(G, V(G)[degree(G) == 0])

  # --- Read the input matrix ---
  if( 'matrix' %in% slotNames(series_object) ){
    matrix_G <- series_object@matrix
  }
  else {
    matrix_G <- read.delim(matrix_path , sep ="\t" ,row.names=1)
  }

  # --- From the peaks in the counts matrix create a GRanges object ---

  peaks <- matrix_G |> rownames() |> as.data.frame() |> dplyr::rename(peaks = 1) |> 
    tidyr::separate('peaks' , c('Chromosome' , 'Start','End') , sep = '-')
  peaks_ranges <- makeGRangesFromDataFrame(peaks)
  peaks_ranges$nodes <- rownames(matrix_G)
  
  network_obj <- new(
    "TCrhoNetNetwork",
    graph = G,
    #threshold = threshold,
    clusters = series_object@communities[[as.character(threshold)]],
    matrix = matrix_G,
    genomicRegions = peaks_ranges,
    modularity = as.data.frame(series_object@modularity[[as.character(threshold)]]),
    metadata = list( n_nodes = vcount(G), n_edges = ecount(G))#source = network_path,
  )

  if (verbose) {
    message(sprintf(
      "✅ Network built: %d nodes, %d edges, transitivity = %.3f",
      vcount(G), ecount(G), series_object@metrics[series_object@metrics$threshold == threshold , 'transitivity']
    ))
  }
  
  return(network_obj)
}