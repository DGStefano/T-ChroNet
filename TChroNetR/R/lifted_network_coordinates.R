#' Lift the coordinates from a genome to another
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
#' @param chain_path Path to a chain file for liftover
#' @param source_build Name of source genome (default = 'hg19')
#' @param target_build Name of target genome (default = 'hg38')
#' @return A TChroNetNetwork object
#' @export
lift_network_coordinates <- function(object, chain_path, source_build = "hg19", target_build = "hg38") {
  if (!inherits(object, "TCrhoNetNetwork")) {
    stop("The 'object' must be a TCrhoNetNetwork instance.")
  }
  if (!file.exists(chain_path)) {
    stop("The UCSC chain file does not exist: ", chain_path)
  }

  # --- 1. Import chain file and perform liftOver ---
  chain <- rtracklayer::import.chain(chain_path)
  gr_hg38_list <- rtracklayer::liftOver(object@genomicRegions, chain)

  # --- 2. Extract original coordinates ---
  orig_coords <- data.frame(
    node = object@genomicRegions$nodes,
    chr_orig = as.character(GenomicRanges::seqnames(object@genomicRegions)),
    start_orig = GenomicRanges::start(object@genomicRegions),
    end_orig = GenomicRanges::end(object@genomicRegions)
  )

  # --- 3. Extract lifted coordinates safely ---
  lifted_coords <- lapply(gr_hg38_list, function(x) {
    if (length(x) == 0) {
      return(data.frame(chr_hg38 = NA, start_hg38 = NA, end_hg38 = NA))
    } else {
      x <- x[1]
      return(data.frame(
        chr_hg38 = as.character(GenomicRanges::seqnames(x)),
        start_hg38 = GenomicRanges::start(x),
        end_hg38 = GenomicRanges::end(x)
      ))
    }
  })

  lifted_df <- do.call(rbind, lifted_coords)

  # --- 4. Handle mismatched lengths ---
  n_nodes <- length(object@genomicRegions)
  if (nrow(lifted_df) != n_nodes) {
    warning("Mismatch between lifted regions and genomicRegions. Adjusting to shortest length.")
    lifted_df <- lifted_df[seq_len(min(nrow(lifted_df), n_nodes)), ]
  }

  # --- 5. Combine original and lifted coordinates ---
  combined_df <- cbind(orig_coords, lifted_df)
  combined_df$lifted_coord <- ifelse(
    is.na(combined_df$chr_hg38),
    NA,
    paste0(combined_df$chr_hg38, "-", combined_df$start_hg38, "-", combined_df$end_hg38)
  )

  # --- 6. Store in object slot ---
  object@lifted_coords <- combined_df

  message(
    "✅ LiftOver complete: ",
    sum(!is.na(combined_df$chr_hg38)), " / ", nrow(combined_df),
    " regions successfully mapped from ", source_build, " → ", target_build, "."
  )

  return(object)
}
