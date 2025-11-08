#' Plot the modularity. If TCrhoNetSeries object is provided the plot represents all the available threshold
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
#' @return A ggplot2 object
#' @export
#' 
plot_modularity <- function(object){
   if (!inherits(object, "TCrhoNetNetwork") & !inherits(object, "TCrhoNetSeries")) {
    stop("The input must be a 'TCrhoNetNetwork' object or 'TCrhoNetSeries'")
  }

  if(inherits(object, "TCrhoNetNetwork")) {
    g <- ggplot(object@modularity , aes (x = resolution , y = modularity , group = 1))+
    geom_line() +
    geom_point()+
    theme_classic()
  }
  else {
    g <- dplyr::bind_rows(object@modularity) |> 
      ggplot(aes (x = resolution , y = modularity , group = threshold , color=threshold))+
      geom_line() +
      geom_point()+
      theme_classic()
  }
  
  return(g)

}