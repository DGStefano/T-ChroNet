#' Plot stacked barplot of annotations per communities
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
#' @param resolutions Set the resolution parameter for plotting trends. If not specified the function the resolution parameter stored in the object (deafault = NULL)
#' @return A ggplot2 object
#' @export
#' 
plot_stacked_annotation <- function(object , resolution = NULL) {
   if (is.null(resolution)){
    resolution = object@resolution
  }

  annotation_df <- genomic_position_stackbar(object , resolution = resolution)

  stacked_data <- annotation_df |> 
    group_by(community_number, annotation) %>%
    summarise(count = n(), .groups = "drop") %>%
    group_by(community_number) %>%
    mutate(proportion = count / sum(count))
  g <- ggplot(stacked_data, aes(fill=annotation, y=proportion, x=community_number)) + 
      geom_bar(position="fill", stat="identity")+
      theme_minimal()
  return(g)
}