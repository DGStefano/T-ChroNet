#' Find the best threshold value according to the transitivity values at increased thresholds
#' 
#' @param object TCrhoNetSeries object
#' @return A TCrhoNetSeries object
#' @export
find_best_th <- function(object , parameter = "relative_lcc") {
  if (!inherits(object, "TCrhoNetSeries")) {
    stop("Input must be a 'TCrhoNetSeries' object.")
  }
 
  if (parameter == "relative_lcc"){
    # Extract metrics
    metrics <- object@metrics
    # Find maximum transitivity
    max_trans <- max(metrics$relative_lcc, na.rm = TRUE)
    # Get indices of rows with max transitivity
    idx <- which(metrics$relative_lcc == max_trans)
    # Pick the highest threshold among tied maxima
    best_th <- max(as.numeric(metrics$threshold[idx]))
  }
  else {
    # Extract metrics
    metrics <- object@metrics
    # Find maximum transitivity
    max_trans <- max(metrics$transitivity, na.rm = TRUE)
    # Get indices of rows with max transitivity
    idx <- which(metrics$transitivity == max_trans)
    # Pick the highest threshold among tied maxima
    best_th <- max(as.numeric(metrics$threshold[idx]))

  }

  # Assign to object
  object@best_th <- best_th

  return(object)
}