# --- ARGUMENT PARSING ---
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript wgcna_bench.R <threads> <input_matrix_path>")
}

n_threads   <- as.numeric(args[1])
input_path  <- args[2]

# --- REMOVE GARBAGE ---
gc(full = TRUE)
# --- FORCE WGNA TO USE FIXED NUMBER OF CORES ---
Sys.setenv(OMP_NUM_THREADS = n_threads)
Sys.setenv(MKL_NUM_THREADS = n_threads)
Sys.setenv(OPENBLAS_NUM_THREADS = n_threads)
Sys.setenv(VECLIB_MAXIMUM_THREADS = n_threads)

# This controls the R-level parallelization
if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    RhpcBLASctl::blas_set_num_threads(n_threads)
    RhpcBLASctl::omp_set_num_threads(n_threads)
}

library(WGCNA)
library(tidyverse)
library(flashClust)

# Enable WGCNA multithreading
if(n_threads > 1) {
    allowWGCNAThreads(nThreads = n_threads)
} else {
    WGCNA::disableWGCNAThreads()
}
# --- DATA LOADING ---
# WGCNA expects genes/peaks in columns, samples in rows
input_mat <- read_tsv(input_path) %>% 
  column_to_rownames(colnames(.)[1]) %>% 
  as.matrix()
input_mat <- t(input_mat)

# --- CORE WGCNA WORKFLOW ---

# 1. Soft Thresholding (Selection of power)
powers <- c(c(1:10), seq(from = 12, to = 20, by = 2))
sft <- pickSoftThreshold(
  input_mat,
  powerVector = powers,
  verbose = 0, # Quiet for benchmarking
  moreNetworkConcepts = TRUE
)

# Use a fixed power if the search fails, or the best one
softPower <- sft$powerEstimate
if(is.na(softPower)) softPower <- 1 

# 2. Calculation of Adjacency and TOM
# This is the heavy computational part
adjacency <- adjacency(input_mat, power = softPower)
TOM <- TOMsimilarity(adjacency, TOMType = "unsigned", verbose = 0)
dissTOM <- 1 - TOM

geneTree = flashClust(as.dist(dissTOM),method="ward");
set.seed(123)
dynamicMods_100 = cutreeDynamic(dendro = geneTree,  method="tree", minClusterSize = 1000)