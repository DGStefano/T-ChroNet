library(tidyverse)
library(ComplexHeatmap)
library(ggplot2)
library(reshape2)

#THP-1
logCPMs_corrected_filtered = read_delim("/mnt/nas-safu/analysis/PhDsdigiove/method_coAcces/data/CellReport/out_mutlicov_matrix/lognorm_edgeR_limma_countsInCellReport_TheiCounts_All.tsv" , delim = "\t") |> column_to_rownames('peaks')

for (file_name in list.files("/mnt/nas-safu/analysis/PhDsdigiove/method_coAcces/data/CellReport/communities_k27/res_1_5/beds/")) {
  comm_file_name = paste("/mnt/nas-safu/analysis/PhDsdigiove/method_coAcces/data/CellReport/communities_k27/res_1_5/beds/" , file_name , sep = "")
  comm_file = read_delim(comm_file_name , col_names = c('chromosome','start','end')) |>
    unite('peaks', c('chromosome', 'start', 'end'), sep = "-") 

  intersected_counts  <- logCPMs_corrected_filtered[comm_file$peaks,]
  intersected_counts.scaled <- t(scale(t(intersected_counts)))


col_fun = circlize::colorRamp2(c(-2,-1,0,1,2), c("blue", "#3f8ef5","white","#e8cc81","#ee2326"))
r.cluster=hclust(stats::dist(intersected_counts.scaled, method = "euclidean"), method="ward.D2")
c.cluster=hclust(stats::dist(t(intersected_counts.scaled), method = "euclidean"), method="ward.D2")

HM = ComplexHeatmap::Heatmap(
  name = "Zscore",
  as.matrix(intersected_counts.scaled),
  # col=col_fun,
  show_column_names = T,
  show_row_names = FALSE,
  cluster_rows = TRUE,
  cluster_columns =FALSE,
  column_order = colnames(intersected_counts.scaled),
  # row_split = 5,
  #column_split = 2,
  na_col = "white",
  row_gap = unit(2, "mm"),
  column_gap = unit(2, "mm"),
  row_title = " ",
  column_title = " ",
  # rect_gp = gpar(col="white", lwd=2),
  column_names_gp = grid::gpar(fontsize = 10),
  row_names_gp = grid::gpar(fontsize = 10),
  # row_dend_reorder = TRUE,
  # show_row_dend = TRUE,
  
)
out_dir = "/mnt/nas-safu/analysis/PhDsdigiove/method_coAcces/data/CellReport/pictures/heatmaps_2_script_final/"
png(paste0(out_dir ,file_name , '.png' , sep ="" ) , height = 5 , width = 9 , units = 'in' , res = 300)
ComplexHeatmap::draw(HM)
dev.off()
}

# LIVER DEVELOPMENT
logCPMs_corrected_filtered = read_delim("/mnt/nas-safu/analysis/PhDsdigiove/method_coAcces/data/liver_devel_mouse_ENCODE/network_files/normalized_samplemean_multicov_all_sites_all_timepoint.tsv" , delim = "\t") |> column_to_rownames('peaks')

for (file_name in list.files("/mnt/nas-safu/analysis/PhDsdigiove/method_coAcces/data/BALL/communities_broad/res_0_9/beds/")) {
  comm_file_name = paste("/mnt/nas-safu/analysis/PhDsdigiove/method_coAcces/data/BALL/communities_broad/res_0_9/beds/" , file_name , sep = "")
  comm_file = read_delim(comm_file_name , col_names = c('chromosome','start','end')) |>
    unite('peaks', c('chromosome', 'start', 'end'), sep = "-") 

  intersected_counts  <- logCPMs_corrected_filtered[comm_file$peaks,]
  intersected_counts.scaled <- t(scale(t(intersected_counts)))


  col_fun = circlize::colorRamp2(c(-2,-1,0,1,2), c("blue", "#3f8ef5","white","#e8cc81","#ee2326"))
  r.cluster=hclust(stats::dist(intersected_counts.scaled, method = "euclidean"), method="ward.D2")
  c.cluster=hclust(stats::dist(t(intersected_counts.scaled), method = "euclidean"), method="ward.D2")

  HM = ComplexHeatmap::Heatmap(
    name = "Zscore",
    as.matrix(intersected_counts.scaled),
    # col=col_fun,
    show_column_names = T,
    show_row_names = FALSE,
    cluster_rows = TRUE,
    cluster_columns =FALSE,
    column_order = colnames(intersected_counts.scaled),
    # row_split = 5,
    #column_split = 2,
    na_col = "white",
    row_gap = unit(2, "mm"),
    column_gap = unit(2, "mm"),
    row_title = " ",
    column_title = " ",
    # rect_gp = gpar(col="white", lwd=2),
    column_names_gp = grid::gpar(fontsize = 10),
    row_names_gp = grid::gpar(fontsize = 10),
    # row_dend_reorder = TRUE,
    # show_row_dend = TRUE,
    
  )
  out_dir = "/mnt/nas-safu/analysis/PhDsdigiove/method_coAcces/data/liver_devel_mouse_ENCODE/pictures/heatmaps_2_script_final/"
  png(paste0(out_dir ,file_name , '.png' , sep ="" ) , height = 5 , width = 9 , units = 'in' , res = 300)
  ComplexHeatmap::draw(HM)
  dev.off()
}

# BCP-ALL
logCPMs_corrected_filtered = read_delim("/mnt/nas-safu/analysis/PhDsdigiove/method_coAcces/data/BALL/allPatientsPatientsNOMERGED_Ball_Multicov_nfcore.tsv" , delim = "\t") |> column_to_rownames('peaks')

for (file_name in list.files("/mnt/nas-safu/analysis/PhDsdigiove/method_coAcces/data/BALL/communities_broad/res_0_9/beds/")) {
  comm_file_name = paste("/mnt/nas-safu/analysis/PhDsdigiove/method_coAcces/data/BALL/communities_broad/res_0_9/beds/" , file_name , sep = "")
  comm_file = read_delim(comm_file_name , col_names = c('chromosome','start','end')) |>
    unite('peaks', c('chromosome', 'start', 'end'), sep = "-") 

  intersected_counts  <- logCPMs_corrected_filtered[comm_file$peaks,]
  intersected_counts.scaled <- t(scale(t(intersected_counts)))
  print(nrow(intersected_counts))

  col_fun = circlize::colorRamp2(c(-2,-1,0,1,2), c("blue", "#3f8ef5","white","#e8cc81","#ee2326"))
  r.cluster=hclust(stats::dist(intersected_counts.scaled, method = "euclidean"), method="ward.D2")
  # c.cluster=hclust(stats::dist(t(intersected_counts.scaled), method = "euclidean"), method="ward.D2")
  c.ha <- ComplexHeatmap::HeatmapAnnotation( CancerStatus = rep(c('Healthy','Primary','Remission','Relapse') , c(6, 11 , 8 , 8) ) ,
                            col = list(CancerStatus = c("Healthy" = "#98C9BF", "Primary" = "#EECA76", "Remission" = "#BC6F61" , 'Relapse' = "#73394F")),
                            na_col = 'white')


  HM = ComplexHeatmap::Heatmap(
    name = "Zscore",
    as.matrix(intersected_counts.scaled),
    col=col_fun,
    show_column_names = FALSE,
    show_row_names = FALSE,
    cluster_rows = r.cluster,
    cluster_columns =FALSE,
    top_annotation = c.ha,
    #column_order = colnames(intersected_counts.scaled),
    # row_split = 5,
    #column_split = 2,
    na_col = "white",
    row_gap = unit(2, "mm"),
    column_gap = unit(2, "mm"),
    row_title = " ",
    column_title = " ",
    # rect_gp = gpar(col="white", lwd=2),
    column_names_gp = grid::gpar(fontsize = 10),
    row_names_gp = grid::gpar(fontsize = 10),
    # row_dend_reorder = TRUE,
    # show_row_dend = TRUE,
    
  )
  out_dir = "/mnt/nas-safu/analysis/PhDsdigiove/method_coAcces/data/BALL/communities_broad/res_0_9/heatmaps_2_script_final/"
  png(paste0(out_dir ,file_name , '.png' , sep ="" ) , height = 5 , width = 9 , units = 'in' , res = 300)
  ComplexHeatmap::draw(HM)
  dev.off()
}

# BCPALL median trends
# Patient groups and their sizes
patient_groups <- c(Healthy = 6, Primary = 11, Remission = 8, Relapse = 8)
group_order <- c("Healthy", "Primary", "Remission", "Relapse")

# Load data
data_matrix_df <- read_delim("/mnt/nas-safu/analysis/PhDsdigiove/method_coAcces/data/BALL/allPatientsPatientsNOMERGED_Ball_Multicov_nfcore.tsv" , delim = "\t") |> column_to_rownames('peaks')


# Generate value_vars (column names of patients)
id_vars <- data_matrix_df$peaks
value_vars <- setdiff(colnames(data_matrix_df), "peaks")

# Create group mapping
group_mapping <- c()
start_idx <- 1
for (group_name in names(patient_groups)) {
  size <- patient_groups[[group_name]]
  group_mapping[value_vars[start_idx:(start_idx + size - 1)]] <- group_name
  start_idx <- start_idx + size
}

# Loop over communities
for (file_name in list.files("/mnt/nas-safu/analysis/PhDsdigiove/method_coAcces/data/BALL/communities_broad/res_0_9/beds/")) {
  comm_file_name = paste("/mnt/nas-safu/analysis/PhDsdigiove/method_coAcces/data/BALL/communities_broad/res_0_9/beds/" , file_name , sep = "")
  comm_file = read_delim(comm_file_name , col_names = c('chromosome','start','end')) |>
    unite('peaks', c('chromosome', 'start', 'end'), sep = "-") 

  # Subset to the current community
  subset_df <- zscored_df[rownames(zscored_df) %in% comm_file$peaks  , ]
  
  # Melt to long format
  melted <- melt(subset_df , variable.name = "original_column", value.name = "Zscore")
  
  # Map group
  melted$group <- group_mapping[as.character(melted$original_column)]
  melted$group <- factor(melted$group, levels = group_order, ordered = TRUE)
  
  # Plot
  p <- ggplot(melted, aes(x = group, y = Zscore)) +
    stat_summary(fun = mean, geom = "line", aes(group = 1), color = "steelblue") +
    ylim(-1, 1) +
    labs(y = "Z-score", x = NULL) +
    theme_classic()

  # Save
  path_to_save <- paste0("/mnt/nas-safu/analysis/PhDsdigiove/method_coAcces/data/BALL/communities_broad/res_0_9/communities_trends_median_final/", file_name,".png")
  ggsave(path_to_save, plot = p, width = 9, height = 5, dpi = 300)
}
