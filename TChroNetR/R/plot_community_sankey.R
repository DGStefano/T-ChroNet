#' Plot the sankey plot of communities at increased resolution for a fixed threshold
#' @import ggsankey
#' @import ggplot2
#' @import dplyr
#' 
#' 
#' @param object TCrhoNetSeries object
#' @param threshold Set the threshold 
#' @return A ggplot2 object
#' @export
plot_community_sankey <- function(object, threshold, output_path = NULL) {
  
  # --- Check input ---
  if (!inherits(object, "TCrhoNetSeries") & !inherits(object, "TCrhoNetNetwork") ) {
    stop("Input must be a 'TCrhoNetSeries' or 'TCrhoNetNetwork' object.")
  }
  if(inherits(object, "TCrhoNetSeries")) {
    if (!threshold %in% names(object@communities)) {
      stop(paste("Threshold", threshold, "not found in object@communities."))
    }
    # --- Extract the data ---
    d <- object@communities[[as.character(threshold)]]
    
    if (!"node" %in% colnames(d)) {
      stop("Community table must contain a 'node' column.")
    }
    
    # Convert to long format for Sankey plotting
    d_long <- d %>% 
      ggsankey::make_long(colnames(d)[colnames(d) != "node"])
    
    # --- Plot ---
    p <- ggplot(d_long, 
                aes(x = x, 
                    next_x = next_x, 
                    node = node, 
                    next_node = next_node,
                    fill = factor(node))) +
      geom_sankey(alpha = 0.8, show.legend = FALSE) +
      scale_fill_discrete(drop = FALSE) +
      theme_bw(base_size = 14) +
      labs(
        title = paste("Sankey Plot of Communities – Threshold", threshold),
        x = "Resolution",
        y = "Community Flow"
      )
    
    # --- Optionally save ---
    if (!is.null(output_path)) {
      ggsave(output_path, p, height = 8, width = 12, units = "in", dpi = 300)
      message("✅ Sankey plot saved to: ", output_path)
    }
  }
  else {
    d <- object@clusters
    # Convert to long format for Sankey plotting
    d_long <- d %>% 
      ggsankey::make_long(colnames(d)[colnames(d) != "node"])
    
    # --- Plot ---
    p <- ggplot(d_long, 
                aes(x = x, 
                    next_x = next_x, 
                    node = node, 
                    next_node = next_node,
                    fill = factor(node))) +
      geom_sankey(alpha = 0.8, show.legend = FALSE) +
      scale_fill_discrete(drop = FALSE) +
      theme_minimal(base_size = 14) +
      theme(axis.text.y = element_blank()) +
      geom_vline(
        xintercept = seq_along(levels(d_long$x)),
        linetype = "dashed",
        color = "grey70" ) +
      labs(
        title = paste("Sankey Plot of Communities – Threshold", threshold),
        x = "Resolution",
        y = "Community Flow"
      )
      
    # --- Optionally save ---
    if (!is.null(output_path)) {
      ggsave(output_path, p, height = 8, width = 12, units = "in", dpi = 300)
      message("✅ Sankey plot saved to: ", output_path)
    }
  }
  
  return(p)
}