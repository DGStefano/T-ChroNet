#' Plot communities trends at a specified resolution
#' 
#' @import ggplot2
#' @importFrom dplyr bind_rows
#' @importFrom tidyr pivot_longer
#' @importFrom scales squish
#' 
#' @param object TChroNetNetwork object
#' @param resolutions Set the resolution parameter for plotting trends. If not specified the function the resolution parameter stored in the object (deafault = NULL)
#' @param custom_ylim list of min and max value for plotted trends (deafult = c(-2,2))
#' @return A ggplot2 object
#' @export

plot_trends_zscore <- function(object,
                               resolution = NULL,
                               custom_ylim = c(-2, 2),
                               time_order = NULL) {
  
  if (!inherits(object, "TCrhoNetNetwork")) {
    stop("The 'object' must be a TCrhoNetNetwork instance.")
  }
  
  if (is.null(resolution)) {
    resolution <- object@resolution
  }
  
  cluster_col <- paste0("clusters_", resolution)
  if (!(cluster_col %in% colnames(object@clusters))) {
    stop(paste0("No clustering results for resolution ", resolution))
  }
  
  clusters_df <- object@clusters[, c("node", cluster_col)]
  colnames(clusters_df) <- c("node", "cluster")
  
  membership_nodes_list <- split(clusters_df$node, clusters_df$cluster)
  
  # ---- compute z-score per peak ----
  count_mat <- object@matrix
  count_df_scored <- t(scale(t(as.matrix(count_mat))))
  
  # ---- long-format table ----
  final_df <- lapply(seq_along(membership_nodes_list), function(i) {
    peaks <- membership_nodes_list[[i]]
    sub <- count_df_scored[rownames(count_df_scored) %in% peaks,, drop = FALSE]

    if (nrow(sub) == 0) return(NULL)

    df <- as.data.frame(sub)
    df$peak <- rownames(df)

    df_long <- tidyr::pivot_longer(
      df, cols = -peak,
      names_to = "time",
      values_to = "zscore"
    )
    df_long$community <- factor(i)
    df_long
  }) |> dplyr::bind_rows()

  if (nrow(final_df) == 0)
    stop("No peaks overlap between matrix and community nodes.")
  
  # ---- ordering of time axis ----
  if (is.null(time_order)) {
    # extract numeric time values from column names
    time_numeric <- suppressWarnings(as.numeric(final_df$time))
    if (all(!is.na(time_numeric))) {
      time_order <- sort(unique(time_numeric))
      final_df$time <- factor(time_numeric, levels = time_order)
    } else {
      final_df$time <- factor(final_df$time, levels = unique(final_df$time))
    }
  } else {
    final_df$time <- factor(final_df$time, levels = time_order)
  }

  # ---- plotting ----  
  g <- ggplot(final_df, aes(x = time, y = zscore)) +
    geom_violin(
      aes(group = time),
      # fill = "gry",
      color = "black",
      alpha = 0.55,
      scale = "width",
      linewidth = 0.3
    ) +
    geom_boxplot(,
  outliers = F,
  width=0.2)+
    stat_summary(
      fun = mean,
      geom = "line",
      aes(group = 1),
      linewidth = 1,
      color = "black"
    ) +
    facet_wrap(~ community, ncol = 1, strip.position = "left") +
    scale_y_continuous(limits = custom_ylim, oob = scales::squish) +
    labs(
      x = "time (min)",
      y = "Z-score (TPM)"
    ) +
    theme_bw(base_size = 12) +
    theme(
      strip.background = element_blank(),
      strip.placement = "outside",
      strip.text = element_text(size = 12, face = "bold"),
      axis.text.x = element_text(angle = 90, hjust = 1),
      legend.position = "none",
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", linewidth = 0.4),
      panel.spacing = unit(1.2, "lines")
    )
  
  return(g)
}
