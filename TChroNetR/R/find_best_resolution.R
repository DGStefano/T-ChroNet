#' Find the best resolution paramenter according to evaluated modularities
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
#' @return An object of class TChroNetNetwork
#' @export
find_best_resolution <- function (object){
  
  if (!inherits(object, "TCrhoNetNetwork")) {
    stop("The 'object' must be a TCrhoNetNetwork instance.")
  }
  
  object@resolution <- object@modularity$resolution[which.max(object@modularity$modularity)]
  

  message(sprintf("Setting best resolution according the modularty at %.2f...", object@resolution))
  return(object)
}