#' Find the best threshold value according to the transitivity values at increased thresholds
#' 
#' @param object TCrhoNetSeries object
#' @return A TCrhoNetSeries object
#' @export
find_best_th <- function(object){
  if (!inherits(object, "TCrhoNetSeries")) {
    stop("Input must be a 'TCrhoNetSeries' object.")
  }
  object@best_th <- as.numeric(object@metrics$threshold[which.max(object@metrics$transitivity)])

  return(object)
}