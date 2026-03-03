#' Interrogate GREAT and retrieve a list of target genes
#'
#' @importFrom rGREAT submitGreatJob getEnrichmentTables getRegionGeneAssociations
#' 
#' @param gr GRanges object
#' @param comm_name Name of the community
#' @param species Name of the genome to be used
#' @return A list of target genes
#' @export
interrogate_GREAT_regions <- function(gr,
                                      comm_name,
                                      species) {
  suppressPackageStartupMessages({
    library(rGREAT)
    library(dplyr)
  })
  
  job <- submitGreatJob(
    gr,
    version = "4.0.4",
    species = species,
    rule = "basalPlusExt",
    adv_upstream = 2.0,
    adv_downstream = 1.0,
    request_interval = 12
  )
  
  tbl <- getEnrichmentTables(job)
  
  target_genes <- getRegionGeneAssociations(job)
  
  return(target_genes)
}