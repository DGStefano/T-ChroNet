#' Convert a TCrhoNetSeries into a TCrhoNetNetwork 
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
  
  if (missing(threshold) & !('best_rh' %in% slotNames(series_object))  ) {
    stop("❌ You must specify a 'threshold' (e.g., 0.7). or run 'find_best_th'")
  }
  if(missing(threshold) & ('best_rh' %in% slotNames(series_object) ) ){
    threshold = series_object@best_rh
  }
  if(missing(matrix_path) ){
    stop("❌ You must specify the matrix path")
  }
  # --- Extract thresholds available ---
  all_thresholds <- as.numeric(names(series_object@communities))
  all_thresholds <- sort(all_thresholds, decreasing = TRUE)
  
  if (!threshold %in% all_thresholds) {
    stop(paste0("❌ Threshold ", threshold, " not found in the series. Available thresholds: ",
                paste(all_thresholds, collapse = ", ")))
  }
  
  # --- Find which edge files are below or equal to the chosen threshold ---
  file_paths <- series_object@edge_files
  file_thrs <- as.numeric(gsub(".*edges_(\\d+\\.\\d+)_.*\\.h5", "\\1", basename(file_paths)))
  
  valid_files <- file_paths[file_thrs >= threshold]
  if (length(valid_files) == 0) stop("❌ No edge files found for this threshold or above.")
  
  valid_files <- valid_files[order(-file_thrs[file_thrs >= threshold])]
  if (verbose) message("📂 Loading ", length(valid_files), " edge files up to threshold ", threshold, "...")
  
  # --- Load and combine edge data ---
  cumulative_edges <- data.frame()
  
  for (f in valid_files) {
    if (verbose) message("🔹 Reading ", basename(f))
    
    e <- tryCatch({
      h5read(f, "/")
    }, error = function(e) {
      warning("⚠️ Could not read ", f, ": ", e$message)
      return(NULL)
    })
    
    if (is.null(e)) next
    e_df <- as.data.frame(e)
    if (ncol(e_df) < 3) next
    colnames(e_df) <- c("id", "origin", "target", "corr")[1:ncol(e_df)]
    
    e_df <- e_df[, c("origin", "target", "corr")]
    e_df <- e_df %>% mutate(corr = as.numeric(corr))
    cumulative_edges <- bind_rows(cumulative_edges, e_df)
  }
  
  # --- Remove duplicate edges ---
  cumulative_edges <- cumulative_edges %>%
    distinct(origin, target, .keep_all = TRUE)
  
  # --- Build graph ---
  if (verbose) message("🧩 Building cumulative network graph...")
  G <- graph_from_data_frame(cumulative_edges, directed = FALSE)
  E(G)$weight <- cumulative_edges$corr
  
  # --- Create TCrhoNetNetwork object ---
  if (verbose) message("🧠 Creating TCrhoNetNetwork object at threshold ", threshold, "...")
  
  # --- Read the input matrix ---
  matrix_G <- read.delim(matrix_path , sep ="\t" ,row.names=1)

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
      vcount(G), ecount(G), series_object@metrics[sereis_network@metrics$threshold == threshold , 'transitivity']
    ))
  }
  
  return(network_obj)
}