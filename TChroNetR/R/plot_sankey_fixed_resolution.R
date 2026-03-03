#' Plot the sankey plot of communities at increased threshod for a fixed resolution
#' 
#' @import ggplot2
#' @importFrom ggsankey make_long geom_sankey
#' @importFrom rlang sym
#' @importFrom dplyr select rename pull
#' @importFrom tibble as_tibble
#' @importFrom methods slotNames
#' 
#' @param object TCrhoNetSeries object
#' @param resolution Fixed resolution to compare
#' @param thresholds List of threshold to compare. If not provided it uses all the sotred thresholds
#' @return A ggplot2 object
#' @export
plot_sankey_fixed_resolution <- function(object,
                                         resolution = NULL,
                                         thresholds = NULL,
                                         output_path = NULL) {

  # --- Check input ---
  if (!inherits(object, "TCrhoNetSeries")) {
    stop("Input must be a 'TCrhoNetSeries' object.")
  }
  
  if (is.null(thresholds)) {
    thresholds <- names(object@communities)
  }
  if(missing(resolution) & ('resolution' %in% slotNames(object) ) ){
    resolution = object@resolution
  }
  
  thresholds <- as.character(sort(as.numeric(thresholds), decreasing = FALSE))
  res_scc_fixed <- NULL
  cnames <- c("Node")
  
  # --- Collect cluster assignments at fixed resolution for each threshold ---
  for (th in thresholds) {
    df <- object@communities[[th]]
    
    if (is.null(df)) {
      warning("No community data found for threshold ", th)
      next
    }
    
    cluster_col <- paste0("clusters_", resolution)
    if (!cluster_col %in% colnames(df)) {
      warning("Resolution ", resolution, " not found for threshold ", th)
      next
    }
    
    new_name <- paste0("TH = ", th)
    cls <- df %>%
      dplyr::select(Node = node, !!sym(cluster_col)) |> 
      dplyr::rename(!!new_name := !!sym(cluster_col)) |> 
      dplyr::pull(!!sym(new_name))
    
    res_scc_fixed <- cbind(res_scc_fixed, cls)
    cnames <- c(cnames, new_name)
  }
  
  # --- Assemble full data frame ---
  res_scc_fixed <- cbind(object@communities[[thresholds[1]]]$node, res_scc_fixed)
  res_scc_fixed <- as.data.frame(res_scc_fixed)
  colnames(res_scc_fixed) <- cnames
  res_scc_fixed <- tibble::as_tibble(res_scc_fixed)
  
  # --- Convert to long format for Sankey plotting ---
  res_scc_fixed_long <- res_scc_fixed |> 
    ggsankey::make_long(colnames(res_scc_fixed)[colnames(res_scc_fixed) != "Node"])
  
  # --- Create Sankey plot ---
  p <- ggplot(res_scc_fixed_long,
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
      xintercept = seq_along(levels(res_scc_fixed_long$x)),
      linetype = "dashed",
      color = "grey70" ) +
    labs(
      title = paste("Sankey Plot of Fixed Resolution =", resolution),
      subtitle = paste("Across thresholds:", paste(thresholds, collapse = ", ")),
      x = "Threshold",
      y = "Community Flow"
    )
  
  # --- Optionally save ---
  if (!is.null(output_path)) {
    ggsave(output_path, p, height = 8, width = 12, units = "in", dpi = 300)
    message("✅ Sankey plot saved to: ", output_path)
  }
  
  return(p)
}