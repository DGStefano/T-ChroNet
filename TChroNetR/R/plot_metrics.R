#' Plot the transitivity at each threshold
#' 
#' @import ggplot2
#' @importFrom rlang sym
#' 
#' @param object TCrhoNetSeries object
#' @return A ggplot2 object
#' @export
plot_metrics <- function(object , metric = "relative_lcc"){
  if (!inherits(object, "TCrhoNetSeries")) {
    stop("The 'object' must be a TCrhoNetSeries instance.")
  }

  if (!(metric %in% colnames(object@metrics)) ){
    stop("The 'metric' is not present in metrics dataframe.")
  }

  metric <- sym(metric)
  g <- ggplot(object@metrics , aes(x = threshold , y = !!metric , group = 1)) +
    geom_point() +
    geom_line()+
    theme_classic()

  return(g)
}