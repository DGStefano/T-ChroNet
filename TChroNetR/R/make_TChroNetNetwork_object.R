#' Create a NetworkAnalysis object from an HDF5 network file
#' 
#' @importFrom igraph make_empty_graph add_vertices add_edges vcount ecount
#' @importFrom arrow read_parquet
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom methods new
#' @importFrom dplyr select
#' 
#' @param edge_files List of parquet files to build the network
#' @param matrix_path Path to the tsv normalized counts matrix
#' @param threshold Desired threshold already selected. If specified, the network is built considering the threshold as the lowest correlation. deafulat = 0
#' @return An object of class NetworkAnalysis
#' @export
make_TChroNetNetwork_obj <- function(edge_files , matrix_path  , threshold = 0) {
  
  if (is.na(threshold)) {
    thresholds <- gsub(".*edges_(\\d+\\.\\d+)_.*\\.parquet", "\\1", basename(edge_files))
    thresholds <- sub("(\\d+\\.\\d).*", "\\1", thresholds)
  }
  else {
    thresholds <- gsub(".*edges_(\\d+\\.\\d+)_.*\\.parquet", "\\1", basename(edge_files))
    thresholds <- sub("(\\d+\\.\\d).*", "\\1", thresholds)
    index <- which(thresholds == as.character(threshold))
    thresholds <- thresholds[index:length(thresholds)]
    edge_files <- edge_files[index:length(edge_files)]
    
  }

  # --- Read the input matrix ---
  matrix_G <- read.delim(matrix_path , sep ="\t" ,row.names=1)

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

    # --- Build graph ---
    G <- igraph::add_edges( G , t(e_df[,c(1,2)]) , attr = list(weight = e_df$corr))
  }
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
    modularity = data.frame(),
    metadata = list( n_nodes = vcount(G), n_edges = ecount(G)),
    threshold = threshold
  )
  
  return(obj)
}
  