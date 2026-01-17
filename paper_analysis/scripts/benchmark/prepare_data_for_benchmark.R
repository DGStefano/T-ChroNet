library(edgeR)
library(dplyr)

# Load your master data
master_counts <- read.table("/home/sdigiove/T-ChroNet/paper_analysis/data/BCP-ALL/consensus_peaks.mLb.clN.featureCounts.txt", header=TRUE , skip = 1 ) |> dplyr::select(-Geneid , -Length , -Strand ) |> 
  unite(sites , c('Chr',"Start","End") , sep ="-") 
master_bed <- master_counts |> dplyr::select(sites) 

master_counts <- master_counts |> column_to_rownames("sites")

prepare_subset <- function(n_nodes) {
  # 1. Select the first N peaks (or random selection)
  selected_peak_ids <- sample(rownames(master_counts) , n_nodes)
  
  # 2. Subset Raw Counts
  raw_subset <- master_counts[selected_peak_ids, ]
  
  # 3. Create Normalized Counts (CPM or TMM)
  # WGCNA and your Python tool often perform better on log-transformed normalized data
  dge <- DGEList(counts = raw_subset)
  dge <- calcNormFactors(dge)
  norm_subset <- cpm(dge, log = TRUE, prior.count = 1)
  
  # 4. Subset BED file (for TCseq)
  bed_subset <- master_bed |> dplyr::filter(sites %in% selected_peak_ids)
  
  # 5. Save files
  save_path <- paste0("/home/sdigiove/T-ChroNet/paper_analysis/data/banchmark/counts/data_", n_nodes)
  dir.create(save_path , showWarnings = F)
  raw_subset |> write.table( paste0(save_path, "/raw_counts.tsv"), sep="\t", quote=F , col.names = TRUE)
  norm_subset |> write.table( paste0(save_path, "/norm_counts.tsv"), sep="\t", quote=F, col.names = TRUE)
  bed_subset |> separate("sites" , c("Chr" , "Start" , "End") , sep = "-") |> write.table(paste0(save_path, "/subset.bed"), sep="\t", quote=F, col.names=T, row.names=F)
}

# Run for all your intervals
intervals <- c(10000, 20000, 30000, 40000 , 50000 , 60000, 70000,80000,90000,100000) 
lapply(intervals, prepare_subset)