#' Plot the transitivity at each threshold
#' @import ggplot2
#' 
#' @param object TCrhoNetSeries object
#' @return A ggplot2 object
#' @export
plot_transitivity <- function(object){
  if (!inherits(object, "TCrhoNetSeries")) {
    stop("The 'object' must be a TCrhoNetSeries instance.")
  }

  g <- ggplot(object@metrics , aes(x = threshold , y = transitivity , group = 1)) +
    geom_point() +
    geom_line()+
    theme_classic()

  return(g)
}