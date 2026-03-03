#' @keywords internal
"_PACKAGE"

.onLoad <- function(libname, pkgname) {
  # Ensures S4 methods work correctly across the package
  if (!requireNamespace("methods", quietly = TRUE)) {
    stop("The 'methods' package is required.")
  }
}

.onUnload <- function(libpath) {
  # Clean up any loaded DLLs if you had C++ code (optional)
  library.dynam.unload("YourPackageName", libpath)
}