#' TCrhoNetNetwork class
#'
#' @slot graph igraph object
#' @slot matrix Matrix of normalized log2 counts used for inferring the netowrk
#' @slot genomicRegions GRanges object of genomic regions
#' @slot lifted_coords Lifted genome coordinates from original genome to raget genome
#' @slot clusters Data frame of cluster memberships for different resolution parameters
#' @slot modularity Data frame of modularity values for each resolution
#' @slot metadata Data frame containing the metadata information of netowrk (number of edges, number of nodes, density)
#' @slot resolution Value of resolution to be used. It stores the last resolution value used or the best resolution value
#' @slot annotations Annotation table from BED overlaps
#' @slot GREAT_targets List of GREAT target genes
#' @export
#' 
setClass(
  "TCrhoNetNetwork",
  slots = list(
    graph = "ANY",
    matrix = 'data.frame',
    genomicRegions = "GRanges",
    lifted_coords = "data.frame",
    clusters = "ANY",
    modularity = "data.frame",
    metadata = "list",
    resolution = 'numeric',
    annotations = 'data.frame',
    GREAT_targets = 'list'
  )
)

#' TCrhoNetTCrhoNetSeriesNetwork class
#'
#' @slot edge_files List of edge hd5 files
#' @slot matrix Matrix of normalized log2 counts used for inferring the netowrk
#' @slot thresholds List of all the threshold stored
#' @slot metrics Data frame containing the metrics for each threshold
#' @slot communities List of communities found for each threshold
#' @slot lifted_coords Lifted genome coordinates from original genome to raget genome
#' @slot modularity Data frame of modularity values for each resolution and each threshold
#' @slot metadata Data frame containing the metadata information of netowrk (number of edges, number of nodes, density)
#' @slot components Data frame of detected components at increased resolution
#' @slot method Method used for investigation of communities
#' @slot best_th Value of best threshold value
#' @slot graph Graph with all the threshold loaded
#' @export
#' 
setClass(
  "TCrhoNetSeries",
  slots = list(
    graph = "ANY",
    matrix = 'data.frame',
    edge_files = "character",
    thresholds = "character",
    metrics = "data.frame",
    communities = "list",
    modularity = 'list',
    components = "data.frame",
    method = "character",
    best_th = "numeric"
  )
)

