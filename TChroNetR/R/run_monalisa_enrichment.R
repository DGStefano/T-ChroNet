#' Annotate regions according to a provided bed annotation file 
#'
#' @import dplyr
#' @import BSgenome
#' @import GenomicRanges
#' @import monaLisa
#' @import tidyr
#' @import BiocParallel
#' 
#' @param object TChroNetNetwork object
#' @param pwms PWMatrixList with motifs for which to calculate enrichments.
#' @param resolution Resolution value to calcualte communities enrichemnt, if NULL the default resolution is used if calculated. default = NULL
#' @param genome specify the genome to use, choose among "hg19", "hg38" , "mm10". default = NULL
#' @param num_cores Number of cores to be used in monaLisa analysis. defult = 1
#' @return A TChroNetNetwork object
#' @export
run_monalisa_enrichment <- function(object , pwms , resolution = NULL ,  genome = NULL  , num_cores = 1) {
  if (!inherits(object, "TCrhoNetNetwork")) {
    stop("The 'object' must be a TCrhoNetNetwork instance.")
  }

  if (is.null(resolution)){
    resolution = object@resolution
  }
  # --- 1. Identify the cluster column name ---
  cluster_col <- paste0("clusters_", as.character(resolution))
  
  if (genome == "hg19") {
    genome_obj <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
    species_id <- 9606 # Human
  } else if (genome == "hg38") {
    genome_obj <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
    species_id <- 9606 # Human
  } else if (genome == "mm10") {
    genome_obj <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
    species_id <- 10090 # Mouse
  } else {
    stop("Unsupported genome. Please use 'hg19', 'hg38', or 'mm10'.")
  }

  monaLisa_enrichments <- list()
  for (i in sort(unique(object@clusters[[cluster_col]]))) {
    community_nodes <- object@clusters[object@clusters[cluster_col] == i, ] |> dplyr::select(node) |> tidyr::separate(node,c('chr','start','end') , sep ="-")

    gr <- GenomicRanges::makeGRangesFromDataFrame(
    community_nodes,
    seqnames.field = colnames(community_nodes)[1],
    start.field = colnames(community_nodes)[2],
    end.field = colnames(community_nodes)[3],
    keep.extra.columns = TRUE
    )

    seqs <- BSgenome::getSeq(genome_obj, gr)

    se_genome <- monaLisa::calcBinnedMotifEnrR(
    seqs = seqs,
    pwmL = pwms,                      
    background = "genome",            
    genome = genome_obj, 
    genome.oversample = 2,
    BPPARAM = BiocParallel::MulticoreParam(num_cores)
    )

    monaLisa_enrichments[[i]] <- se_genome
  }
  
  return(monaLisa_enrichments)

}