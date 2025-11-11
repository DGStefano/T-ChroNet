#' @keywords internal
"_PACKAGE"

.onLoad <- function(libname, pkgname) {
  # Optional: Automatically load imported packages
  pkgs <- c("hdf5r","igraph","leidenAlg","GenomicRanges","rtracklayer","liftOver","leidenAlg","dplyr")
  invisible(lapply(pkgs, function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(paste("Package", pkg, "is required but not installed."))
    }
    library(pkg, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)
  }))
}