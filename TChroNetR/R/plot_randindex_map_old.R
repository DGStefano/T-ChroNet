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
plot_randindex_map <- function(object, output_path = NULL , color_limits = c(0.0, 1)) {
  if (!inherits(object, "TCrhoNetSeries")) {
    stop("Input must be a 'TCrhoNetSeries' object.")
  }
  # if (!threshold %in% names(object@communities)) {
  #  stop("Threshold not found in object@communities.")
  # }

  library(fossil)
  library(tidyverse)
  library(grid)

  # --- Extract cluster dataframe for this threshold ---
  # comm_df <- object@communities[[as.character(threshold)]]

  # Remove 'node' column if present
  # if ("node" %in% colnames(comm_df)) {
  #  comm_df <- comm_df %>% tibble::column_to_rownames("node")
  # }

  # --- Prepare a Rand Index matrix structure ---
  n_res <- ncol(object@communities[[1]][2:length(object@communities[[1]])])
  n_th <- length(names(object@communities))
  RandIndex_matrix <- matrix(list(), nrow = n_res, ncol = n_th)
  colnames(RandIndex_matrix) <- names(object@communities)
  rownames(RandIndex_matrix) <- colnames(object@communities[[1]][2:length(object@communities[[1]])])

  # --- Fill RandIndex_matrix with cluster assignment lists ---
  for (i in seq_len(n_res)) {
    for (j in seq_len(n_th)) {
      RandIndex_matrix[[i, j]] <- list(object@communities[[names(object@communities)[j]]][,i+1])
    }
  }

  # --- Function to get adjacent pairs ---
  get_adj_pairs <- function(mat) {
    nr <- nrow(mat)
    nc <- ncol(mat)
    pairs <- data.frame()

    for (i in 1:nr) {
      for (j in 1:nc) {
        if (i < nr) {
          pairs <- rbind(pairs, data.frame(
            col = colnames(mat)[j],
            row = rownames(mat)[i],
            col_comp = colnames(mat)[j],
            row_comp = rownames(mat)[i + 1],
            rand_index = fossil::adj.rand.index(mat[[i, j]][[1]], mat[[i + 1, j]][[1]])
          ))
        }
        if (j < nc) {
          pairs <- rbind(pairs, data.frame(
            col = colnames(mat)[j],
            row = rownames(mat)[i],
            col_comp = colnames(mat)[j + 1],
            row_comp = rownames(mat)[i],
            rand_index = fossil::adj.rand.index(mat[[i, j]][[1]], mat[[i, j + 1]][[1]])
          ))
        }
        if (i < nr & j < nc) {
          pairs <- rbind(pairs, data.frame(
            col = colnames(mat)[j],
            row = rownames(mat)[i],
            col_comp = colnames(mat)[j + 1],
            row_comp = rownames(mat)[i + 1],
            rand_index = fossil::adj.rand.index(mat[[i, j]][[1]], mat[[i + 1, j + 1]][[1]])
          ))
        }
      }
    }
    pairs
  }

  # --- Compute pairwise adj Rand indices ---
  pairs <- get_adj_pairs(RandIndex_matrix)
  print("plotting")
  # --- Plot ---
  p <- ggplot() +
    geom_segment(
      data = pairs,
      aes(
        x = col, y = row,
        xend = col_comp, yend = row_comp,
        color = rand_index
      ),
      arrow = arrow(length = unit(0.2, "inches")),
      alpha = 0.8,
      size = 1
    ) +
    scale_color_viridis_c(option = "C" , limits = color_limits) +
    theme_minimal(base_size = 14) +
    coord_equal() +
    labs(
      title = paste("Adjusted Rand Index Map"),
      color = "Adj. Rand Index",
      x = "Resolution", y = "Resolution"
    )

  if (!is.null(output_path)) {
    ggsave(output_path, p, height = 8, width = 10, units = "in", dpi = 300)
    message("✅ Plot saved to: ", output_path)
  }

  return(p)
}