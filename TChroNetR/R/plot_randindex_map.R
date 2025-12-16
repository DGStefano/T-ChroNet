#' Plot the random index map in which communities are compared among neighbors at increased threshold and resolution.
#' 
#' @import fossil
#' @import dplyr
#' @import grid
#' 
#' @param object TCrhoNetSeries object
#' @param color_limits List of min and max limits for color map (default = c(0.0,1))
#' @return A ggplot2 object
#' @export
#' 
plot_randindex_map <- function(object,
                               output_path = NULL,
                               color_limits = c(0, 1)) {

  if (!inherits(object, "TCrhoNetSeries")) {
    stop("Input must be a 'TCrhoNetSeries' object.")
  }

  comms <- object@communities
  th_names <- names(comms)

  res_names <- colnames(comms[[1]])[-1]
  n_res <- length(res_names)
  n_th  <- length(th_names)

  # Preallocate maximum possible number of edges
  n_edges <- (n_res - 1) * n_th +        # vertical
             n_res * (n_th - 1) +        # horizontal
             (n_res - 1) * (n_th - 1)    # diagonal

  pairs <- data.frame(
    col = character(n_edges),
    row = character(n_edges),
    col_comp = character(n_edges),
    row_comp = character(n_edges),
    rand_index = numeric(n_edges),
    stringsAsFactors = FALSE
  )

  k <- 1

  for (j in seq_len(n_th)) {
    comm_j <- comms[[j]]

    for (i in seq_len(n_res)) {

      v <- comm_j[, i + 1]

      if (i < n_res) {
        pairs[k, ] <- list(
          th_names[j], res_names[i],
          th_names[j], res_names[i + 1],
          fossil::adj.rand.index(v, comm_j[, i + 2])
        )
        k <- k + 1
      }

      if (j < n_th) {
        comm_j1 <- comms[[j + 1]]
        pairs[k, ] <- list(
          th_names[j], res_names[i],
          th_names[j + 1], res_names[i],
          fossil::adj.rand.index(v, comm_j1[, i + 1])
        )
        k <- k + 1
      }

      if (i < n_res && j < n_th) {
        pairs[k, ] <- list(
          th_names[j], res_names[i],
          th_names[j + 1], res_names[i + 1],
          fossil::adj.rand.index(v, comms[[j + 1]][, i + 2])
        )
        k <- k + 1
      }
    }
  }

  pairs <- pairs[seq_len(k - 1), ]
  print("plotting")
  p <- ggplot2::ggplot(pairs) +
    ggplot2::geom_segment(
      ggplot2::aes(
        x = col, y = row,
        xend = col_comp, yend = row_comp,
        color = rand_index
      ),
      arrow = grid::arrow(length = grid::unit(0.2, "inches")),
      alpha = 0.8,
      linewidth = 1
    ) +
    ggplot2::scale_color_viridis_c(limits = color_limits) +
    ggplot2::coord_equal() +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::labs(
      title = "Adjusted Rand Index Map",
      color = "Adj. Rand Index",
      x = "Thresholds",
      y = "Resolution"
    )

  if (!is.null(output_path)) {
    ggplot2::ggsave(output_path, p, height = 8, width = 10, dpi = 300)
  }

  p
}