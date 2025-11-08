#' Plot communities trends at a specified resolution
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
#' @param custom_ylim list of min and max value for plotted trends (deafult = c(-2,2))
#' @return A ggplot2 object
#' @export

plot_trends_zscore <- function(object, resolution = NULL , custom_ylim = c(-2, 2)) {
  if (!inherits(object, "TCrhoNetNetwork")) {
    stop("The 'object' must be a TCrhoNetNetwork instance.")
  }
  
  if (is.null(resolution)){
    resolution = object@resolution
  }

  # --- 1. Identify the cluster column name ---
  cluster_col <- paste0("clusters_", as.character(resolution))
  
  if (!(cluster_col %in% colnames(object@clusters))) {
    stop(paste0("No clustering results found for resolution ", resolution))
  }

  
  # --- 2. Build community list from clusters ---
  clusters_df <- object@clusters[, c("node", cluster_col), drop = FALSE]
  colnames(clusters_df) <- c("node", "cluster")
  
  membership_nodes_list <- split(clusters_df$node, clusters_df$cluster)
  
  # --- 3. Compute Z-score per row ---
  count_df_scored <- t(scale(t(as.matrix(object@matrix))))
  if (any(is.na(count_df_scored))) {
    warning("NA values produced during z-score normalization (check for zero-variance rows).")
  }
  
  # --- 4. Build long-format data for ggplot ---
  final_df <- data.frame()
  
  for (i in seq_along(membership_nodes_list)) {
    peaks <- membership_nodes_list[[i]]
    sub_df <- count_df_scored[rownames(count_df_scored) %in% peaks, , drop = FALSE]
    
    if (nrow(sub_df) == 0) next
    
    sub_df <- as.data.frame(sub_df)
    sub_df$peaks <- rownames(sub_df)
    sub_df <- tidyr::pivot_longer(sub_df, cols = -peaks, names_to = "variable", values_to = "value")
    sub_df$community <- as.factor(i)
    print(sub_df)
    final_df <- rbind(final_df, sub_df)
  }
  
  if (nrow(final_df) == 0) {
    stop("No overlapping peaks between data_matrix rownames and network node names.")
  }
  
  # --- 5. Plot violin + mean line trends ---
  g <- ggplot(final_df, aes(x = variable, y = value, group = community)) +
    geom_violin( aes(x = variable, y = value, group = variable) , scale = "width", fill = 'grey' , alpha = 0.3) +
    geom_line(stat = "summary", fun = mean, color = "black", linewidth = 1) +
    facet_wrap(~ community, ncol = 1, scales = "free_y") +
    scale_y_continuous(limits = custom_ylim, oob = scales::squish) +
    labs(
      x = "Variable",
      y = "Accessibility Z-score",
      title = paste0("Community Accessibility Trends (resolution = ", resolution, ")")
    ) +
    theme_bw(base_size = 12) +
    theme(
      strip.text = element_text(size = 12, face = "bold"),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      legend.position = "none",
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black"),
      plot.title = element_text(hjust = 0.5)
    )
  
  return(g)
}