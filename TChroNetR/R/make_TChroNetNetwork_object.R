#' Create a NetworkAnalysis object from an HDF5 network file
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
#' @param network_path Path to the HDF5 file
#' @param matrix_path Path to the tsv normalized counts matrix
#' @return An object of class NetworkAnalysis
#' @export
make_TChroNetNetwork_obj <- function(network_path , matrix_path ) {
  
  # --- Read the HDF5 file ---
  e <- rhdf5::h5read(network_path, "/")
  e_df <- as.data.frame(e)
  
  # --- Adjust columns (assuming 2nd to 4th are origin, target, corr) ---
  e_df <- e_df[, c(2:ncol(e_df))]
  colnames(e_df) <- c("origin", "target", "corr")
  
  # --- Build the graph ---
  G <- graph_from_data_frame(e_df, directed = FALSE)
  if ("corr" %in% colnames(e_df)) {
    E(G)$corr <- e_df$corr
  }

  # --- Read the input matrix ---
  matrix_G <- read.delim(matrix_path , sep ="\t" ,row.names=1)

  # --- From the peaks in the counts matrix create a GRanges object ---

  peaks <- matrix_G |> rownames() |> as.data.frame() |> dplyr::rename(peaks = 1) |> 
    tidyr::separate('peaks' , c('Chromosome' , 'Start','End') , sep = '-')
  peaks_ranges <- makeGRangesFromDataFrame(peaks)
  peaks_ranges$nodes <- rownames(matrix_G)
  # --- Initialize empty analysis slots ---
  obj <- new(
    "TCrhoNetNetwork",
    graph = G,
    matrix = matrix_G,
    genomicRegions = peaks_ranges,
    clusters =data.frame(),
    modularity = list(),
    metadata = list( n_nodes = vcount(G), n_edges = ecount(G))#source = network_path,
  )
  
  return(obj)
}
