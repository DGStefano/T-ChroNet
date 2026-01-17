#!/home/user/miniconda3/envs/R/bin/Rscript
library(TChroNetR)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
input_path <- args[1]
edge_files <- args[2]

# --- REMOVE GARBAGE ---
gc(full = TRUE)

# List all graph edge files (thresholded networks)
edge_files_list <- list.files(edge_files, full.names = TRUE)

# Path to the normalized count matrix
matrix_path <- input_path

# Build the TChroNet series object
obj_series <- build_TChroNetSeries_object(
  edge_files_list,
  matrix_path,
  run_cd = FALSE,
  transitivity = FALSE
)

obj_series@best_th <- find_best_th(obj_series)

obj <- build_TCrhoNetNetwork_from_series(obj_series )

obj <- compute_membership_nodes(obj,
                                resolutions = seq(0.5, 1.5, 0.1)
                                )
obj <- find_best_resolution(obj)
