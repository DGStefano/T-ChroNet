#' Plot the density at each threshold
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
#' @return A ggplot2 object
#' @export
plot_density <- function(object){
  if (!inherits(object, "TCrhoNetSeries")) {
    stop("The 'object' must be a TCrhoNetSeries instance.")
  }

  g <- ggplot(object@metrics , aes(x = threshold , y = density , group = 1)) +
    geom_point() +
    geom_line()+
    theme_classic()

  return(g)
}