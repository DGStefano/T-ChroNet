library(clustree)
library(tidyverse)
library(readr)
library(tidyverse)
library(ComplexHeatmap)
library(clusterProfiler)
library(msigdbr)
library(TChroNetR)

plot_heatmap_counts <- function(object , resolution = NULL) {
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
  cluster_df <- object@clusters[c("node" , cluster_col)] 
  ordered_rows <- cluster_df[order(cluster_df[[cluster_col]]),]
  ordered_matrix <- object@matrix[ordered_rows[['node']],]
  
  ha_row <- ComplexHeatmap::rowAnnotation( community = as.character(ordered_rows[[cluster_col]])  )
  
  h <- ComplexHeatmap::Heatmap(t(scale(t(ordered_matrix))) ,
  right_annotation = ha_row,
  cluster_rows = F,
  show_row_names = F,
  show_column_dend = F,
  cluster_columns = F
  )
  
  return(h)
}

plot_cistrom <- function(
    base_folder,
    tf_name_file,
    top_n = 10
){

  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(tibble)

  # -------------------------
  # 1. List community folders
  # -------------------------
  folders_to_compare <- list.files(base_folder)

  # final df
  final_homer_tfs <- NULL

  # -------------------------
  # 2. Read each community file
  # -------------------------
  for (community_num in seq_along(folders_to_compare)) {

    community <- folders_to_compare[[community_num]]
    to_read <- file.path(base_folder, community)

    homer_tfs <- read_delim(to_read, delim = ",")

    homer_tfs <- homer_tfs |>
      group_by(Factor) |>
      summarise(max_value = max(GIGGLE_score)) |>
      dplyr::select(Factor, max_value)

    colnames(homer_tfs) <- c("Factor", community_num)

    if (community_num == 1) {
      final_homer_tfs <- homer_tfs
    } else {
      final_homer_tfs <- merge(final_homer_tfs, homer_tfs, by = "Factor", all = TRUE)
    }
  }

  # convert to rownames
  final_homer_tfs <- final_homer_tfs |> column_to_rownames("Factor")

  # -------------------------
  # 3. Load TF name list
  # -------------------------
  tfs_name <- read_delim(tf_name_file, delim = "\t", col_names = FALSE)
  tfs_name <- tfs_name |> mutate(across(where(is.character), toupper))

  # keep only real TFs
  real_tfs <- rownames(final_homer_tfs)[rownames(final_homer_tfs) %in% tfs_name$X1]
  final_homer_tfs_only_tfs <- final_homer_tfs[real_tfs, , drop = FALSE]

  # -------------------------
  # 4. Select top N TFs per community
  # -------------------------
  saved_tfs <- c()

  for (column_name in colnames(final_homer_tfs_only_tfs)) {
    col_df <- final_homer_tfs_only_tfs[!is.na(final_homer_tfs_only_tfs[[column_name]]), column_name, drop = FALSE]
    col_df <- col_df[order(-col_df[[column_name]]),, drop = FALSE]
    top_values <- head(col_df, top_n)
    saved_tfs <- c(saved_tfs, rownames(top_values))
  }

  saved_tfs <- unique(saved_tfs)

  # matrix for ordering
  cluster_map <- final_homer_tfs_only_tfs[saved_tfs, ]
  cluster_map[is.na(cluster_map)] <- 0

  r.cluster <- hclust(dist(cluster_map, method = "euclidean"), method = "ward.D2")

  # -------------------------
  # 5. Long-format for plotting
  # -------------------------
  final_homer_scatterplot <- final_homer_tfs_only_tfs[saved_tfs, ] |>
    rownames_to_column("Factor") |>
    pivot_longer(cols = -Factor, names_to = "Community", values_to = "ComboScore")

  final_homer_scatterplot$Factor <- factor(final_homer_scatterplot$Factor,
                                           levels = r.cluster$labels[r.cluster$order])

  final_homer_scatterplot$Community <- as.numeric(gsub("\\D", "", final_homer_scatterplot$Community))

  final_homer_scatterplot$Community <- factor(
    final_homer_scatterplot$Community,
    levels = sort(unique(final_homer_scatterplot$Community))
)
  # -------------------------
  # 6. Plot
  # -------------------------
  p <- ggplot(final_homer_scatterplot, aes(x = Community, y = Factor)) +
    geom_point(aes(size = ComboScore),
               shape = 21, stroke = 0.6, color = "black") +
    geom_point(aes(size = ComboScore),
               shape = 21, fill = "black", stroke = 0.6, alpha = 0.4) +
    scale_x_discrete(drop=F)+
    theme_classic() +
    ylab("") +
    theme(
      legend.position = "top",
      legend.text = element_text(size = 10)
    )

  return(p)
}

# edge_files_ball <- list.files("/mnt/nas-safu02/sdigiove_workspace/check_th_TCHRONET/new_tchroent_ties/ball/th/" , full.names = T)
# matrix_path <- "~/T-ChroNet/paper_analysis/data/BCP-ALL/allPatientsPatientsNOMERGED_Ball_Multicov_nfcore.tsv"
# series_ball <- build_TChroNetSeries_object(edge_files_ball,matrix_path,method = "Leiden",resolutions_list = seq(0.5, 1.5, 0.1),run_cd = T,seed = 123, transitivity = FALSE)
# write_rds(series_ball , "/mnt/nas-safu02/sdigiove_workspace/check_th_TCHRONET/new_tchroent_ties/ball/series_obj_ball.rds")

series_ball <- read_rds("/mnt/nas-safu02/sdigiove_workspace/check_th_TCHRONET/new_tchroent_ties/ball/series_obj_ball.rds")


pdf("/mnt/nas-safu02/sdigiove_workspace/check_th_TCHRONET/new_tchroent_ties/ball/pictures/rand_index_map.pdf"  ,
  height = 9, 
  width = 5,
  family = "ArialMT",
  useDingbats = FALSE)
plot_randindex_map(series_ball )
dev.off()

pdf("/mnt/nas-safu02/sdigiove_workspace/check_th_TCHRONET/new_tchroent_ties/ball/pictures/relative_lcc.pdf"  ,
  height = 5, 
  width = 9,
  family = "ArialMT",
  useDingbats = FALSE)
plot_metrics(series_ball , "relative_lcc") +
  geom_point( size = 3)+
  xlab("Threshold")+
  ylab("Relative LCC")+
  theme(text = element_text(size = 15))
dev.off()


series_ball@best_th <- find_best_th(series_ball) 

ball_obj <- build_TCrhoNetNetwork_from_series(series_ball)
# write_rds(ball_obj , "/mnt/nas-safu02/sdigiove_workspace/check_th_TCHRONET/new_tchroent_ties/ball/obj_ball.rds")
# ball_obj <- read_rds("/mnt/nas-safu02/sdigiove_workspace/check_th_TCHRONET/new_tchroent_ties/ball/obj_ball.rds")
# colnames(ball_obj@matrix) <- c("H6","H9","H4","H7","H5","H8","P14","P12","P10","P13","P11","P16","P7","P19","P17","P27","P26","REM21","REM22","REM18","REM26","REM20","REM25","REM23","REM24","REL31","REL28","REL29","REL33","REL32","REL27","REL30","REL14")
# ball_obj@clusters <- data.frame()

pdf("/mnt/nas-safu02/sdigiove_workspace/check_th_TCHRONET/new_tchroent_ties/ball/pictures/community_sankey_th_05.pdf"  ,
  height = 5, 
  width = 9,
  family = "ArialMT",
  useDingbats = FALSE)
plot_community_sankey(series_ball , threshold = 0.5)
dev.off()

pdf("/mnt/nas-safu02/sdigiove_workspace/check_th_TCHRONET/new_tchroent_ties/ball/pictures/modularity.pdf"  ,
  height = 5, 
  width = 9,
  family = "ArialMT",
  useDingbats = FALSE)
plot_modularity(ball_obj)
dev.off()


ball_obj <- find_best_resolution(ball_obj)

pdf("/mnt/nas-safu02/sdigiove_workspace/check_th_TCHRONET/new_tchroent_ties/ball/pictures/community_size.pdf"  ,
  height = 5, 
  width = 9,
  family = "ArialMT",
  useDingbats = FALSE)
ball_obj@clusters |> group_by(clusters_1) |> 
  summarise(peaks = n()) |> 
  mutate(clusters_1 = as.factor(clusters_1)) |> 
ggplot(aes(x = clusters_1 , y = peaks)) +
  geom_bar(stat = 'identity' , fill = "black") +
  theme_classic() + 
  coord_flip()+
  xlab("")+
  ylab("Size")+
  theme(text = element_text(size = 15))
dev.off()


pdf("/mnt/nas-safu02/sdigiove_workspace/check_th_TCHRONET/new_tchroent_ties/ball/sankey_plot_best_res.pdf"  ,
  height = 5, 
  width = 9,
  family = "ArialMT",
  useDingbats = FALSE)
plot_sankey_fixed_resolution(series_ball , resolution = 1.0)
dev.off()


ball_obj <- lift_network_coordinates(ball_obj , chain_path = "/mnt/nas-safu01/analysis/scripts/ScriptSdigiove/RegNetATACProject/T-ChroNet/toy_data/data/hg19ToHg38.over.chain")
ball_obj <- annotate_regions_from_bed(ball_obj , bed_path = "/mnt/nas-safu01/analysis/scripts/ScriptSdigiove/RegNetATACProject/T-ChroNet/toy_data/data/GRCh38-cCREs.bed" , genome = 'hg38')
pdf("/mnt/nas-safu02/sdigiove_workspace/check_th_TCHRONET/new_tchroent_ties/ball/stacked.pdf"  ,
  height = 5, 
  width = 9,
  family = "ArialMT",
  useDingbats = FALSE)
plot_stacked_annotation(ball_obj , resolution = as.numeric(ball_obj@resolution)) +
  coord_flip() +
  theme(legend.position = "top" , text = element_text(size = 15)) +
  ylab("") +
  xlab("")
dev.off()


ball_obj <- run_GREAT_analysis(ball_obj, resolution = 1.0 , genome = "hg19") #as.numeric(ball_obj@resolution)

pdf("/mnt/nas-safu02/sdigiove_workspace/check_th_TCHRONET/new_tchroent_ties/ball/pictures/boxplot_trends.pdf"  ,
  height = 10, 
  width = 7,
  family = "ArialMT",
  useDingbats = FALSE)
plot_trends_zscore(ball_obj)+
  theme(text = element_text(size = 15)) +
  scale_x_discrete(
    labels = c(
      "patient_6_Healthy_ATACseq_REP1.mLb.clN.sorted.bam" = "H6",
 "patient_9_Healthy_ATACseq_REP1.mLb.clN.sorted.bam" = "H9",
 "patient_4_Healthy_ATACseq_REP1.mLb.clN.sorted.bam"    = "H4",
 "patient_7_Healthy_ATACseq_REP1.mLb.clN.sorted.bam" = "H7",
 "patient_5_Healthy_ATACseq_REP1.mLb.clN.sorted.bam" = "H5",
 "patient_8_Healthy_ATACseq_REP1.mLb.clN.sorted.bam" = "H8",
 "patient_14_Primary_ATACseq_REP1.mLb.clN.sorted.bam" = "P14",
 "patient_12_Primary_ATACseq_REP1.mLb.clN.sorted.bam" = "P12",
 "patient_10_Primary_ATACseq_REP1.mLb.clN.sorted.bam" = "P10",
 "patient_13_Primary_ATACseq_REP1.mLb.clN.sorted.bam" = "P13",
 "patient_11_Primary_ATACseq_REP1.mLb.clN.sorted.bam" = "P11",
 "patient_16_Primary_ATACseq_REP1.mLb.clN.sorted.bam" = "P16",
 "patient_7_Primary_ATACseq_REP1.mLb.clN.sorted.bam" = "P7",
 "patient_19_Primary_ATACseq_REP1.mLb.clN.sorted.bam" = "P19",
 "patient_17_Primary_ATACseq_REP1.mLb.clN.sorted.bam" = "P17",
 "patient_27_Primary_ATACseq_REP1.mLb.clN.sorted.bam" = "P27",
 "patient_26_Primary_ATACseq_REP1.mLb.clN.sorted.bam" = "P26",
 "patient_21_Remission_ATACseq_REP1.mLb.clN.sorted.bam" = "REM21",
 "patient_22_Remission_ATACseq_REP1.mLb.clN.sorted.bam" = "REM22",
 "patient_18_Remission_ATACseq_REP1.mLb.clN.sorted.bam" = "REM18",
 "patient_26_Remission_ATACseq_REP1.mLb.clN.sorted.bam" = "REM26",
 "patient_20_Remission_ATACseq_REP1.mLb.clN.sorted.bam" = "REM20",
 "patient_25_Remission_ATACseq_REP1.mLb.clN.sorted.bam" = "REM25",
 "patient_23_Remission_ATACseq_REP1.mLb.clN.sorted.bam" = "REM23",
 "patient_24_Remission_ATACseq_REP1.mLb.clN.sorted.bam" = "REM24",
 "patient_31_Relapse_ATACseq_REP1.mLb.clN.sorted.bam" = "REL31",
 "patient_28_Relapse_ATACseq_REP1.mLb.clN.sorted.bam" = "REL28",
 "patient_29_Relapse_ATACseq_REP1.mLb.clN.sorted.bam" = "REL29",
 "patient_33_Relapse_ATACseq_REP1.mLb.clN.sorted.bam" = "REL33",
 "patient_32_Relapse_ATACseq_REP1.mLb.clN.sorted.bam" = "REL32",
 "patient_27_Relapse_ATACseq_REP1.mLb.clN.sorted.bam" = "REL27",
 "patient_30_Relapse_ATACseq_REP1.mLb.clN.sorted.bam" = "REL30",
 "patient_14_Relapse_ATACseq_REP1.mLb.clN.sorted.bam"  = "REL14"
    ))
dev.off()

pdf("/mnt/nas-safu02/sdigiove_workspace/check_th_TCHRONET/new_tchroent_ties/ball/pictures/lineplot_trends.pdf"  ,
  height = 10, 
  width = 7,
  family = "ArialMT",
  useDingbats = FALSE)
plot_line_trends_zscore(ball_obj)+
    theme(text = element_text(size = 15)) +
  scale_x_discrete(
    labels = c(
      "patient_6_Healthy_ATACseq_REP1.mLb.clN.sorted.bam" = "H6",
 "patient_9_Healthy_ATACseq_REP1.mLb.clN.sorted.bam" = "H9",
 "patient_4_Healthy_ATACseq_REP1.mLb.clN.sorted.bam"    = "H4",
 "patient_7_Healthy_ATACseq_REP1.mLb.clN.sorted.bam" = "H7",
 "patient_5_Healthy_ATACseq_REP1.mLb.clN.sorted.bam" = "H5",
 "patient_8_Healthy_ATACseq_REP1.mLb.clN.sorted.bam" = "H8",
 "patient_14_Primary_ATACseq_REP1.mLb.clN.sorted.bam" = "P14",
 "patient_12_Primary_ATACseq_REP1.mLb.clN.sorted.bam" = "P12",
 "patient_10_Primary_ATACseq_REP1.mLb.clN.sorted.bam" = "P10",
 "patient_13_Primary_ATACseq_REP1.mLb.clN.sorted.bam" = "P13",
 "patient_11_Primary_ATACseq_REP1.mLb.clN.sorted.bam" = "P11",
 "patient_16_Primary_ATACseq_REP1.mLb.clN.sorted.bam" = "P16",
 "patient_7_Primary_ATACseq_REP1.mLb.clN.sorted.bam" = "P7",
 "patient_19_Primary_ATACseq_REP1.mLb.clN.sorted.bam" = "P19",
 "patient_17_Primary_ATACseq_REP1.mLb.clN.sorted.bam" = "P17",
 "patient_27_Primary_ATACseq_REP1.mLb.clN.sorted.bam" = "P27",
 "patient_26_Primary_ATACseq_REP1.mLb.clN.sorted.bam" = "P26",
 "patient_21_Remission_ATACseq_REP1.mLb.clN.sorted.bam" = "REM21",
 "patient_22_Remission_ATACseq_REP1.mLb.clN.sorted.bam" = "REM22",
 "patient_18_Remission_ATACseq_REP1.mLb.clN.sorted.bam" = "REM18",
 "patient_26_Remission_ATACseq_REP1.mLb.clN.sorted.bam" = "REM26",
 "patient_20_Remission_ATACseq_REP1.mLb.clN.sorted.bam" = "REM20",
 "patient_25_Remission_ATACseq_REP1.mLb.clN.sorted.bam" = "REM25",
 "patient_23_Remission_ATACseq_REP1.mLb.clN.sorted.bam" = "REM23",
 "patient_24_Remission_ATACseq_REP1.mLb.clN.sorted.bam" = "REM24",
 "patient_31_Relapse_ATACseq_REP1.mLb.clN.sorted.bam" = "REL31",
 "patient_28_Relapse_ATACseq_REP1.mLb.clN.sorted.bam" = "REL28",
 "patient_29_Relapse_ATACseq_REP1.mLb.clN.sorted.bam" = "REL29",
 "patient_33_Relapse_ATACseq_REP1.mLb.clN.sorted.bam" = "REL33",
 "patient_32_Relapse_ATACseq_REP1.mLb.clN.sorted.bam" = "REL32",
 "patient_27_Relapse_ATACseq_REP1.mLb.clN.sorted.bam" = "REL27",
 "patient_30_Relapse_ATACseq_REP1.mLb.clN.sorted.bam" = "REL30",
 "patient_14_Relapse_ATACseq_REP1.mLb.clN.sorted.bam"  = "REL14"
    ))
dev.off()

patient_groups <- c(
  Healthy = 6,
  Primary = 11,
  Remission = 8,
  Relapse = 8
)

order <- c("Healthy", "Primary", "Remission", "Relapse")
value_vars <- colnames(ball_obj@matrix)
group_mapping <- rep(names(patient_groups), times = patient_groups)
names(group_mapping) <- value_vars

facet_colors <- c(
  Healthy   = "#66c5cc",
  Primary   = "#ff6767",
  Remission = "#ffdb84",
  Relapse   = "#ff951f"
)
library(stats)
scaled_matrix <- t(scale(t(ball_obj@matrix))) |> as.data.frame()
all_melted <- lapply(seq_along(list_communities_11), function(x) {
  df <- scaled_matrix[list_communities_11[[x]], ] |> 
    rownames_to_column("peaks") |>
    pivot_longer(
      cols = all_of(value_vars),
      names_to = "original_column",
      values_to = "expression"
    ) |>
    mutate(
      group = factor(group_mapping[original_column], levels = order),
      community = paste0("Community_", x)
    )
  return(df)
}) |> bind_rows()

median_df <- all_melted |> 
  group_by(community, group) |> 
  summarise(median_expr = median(expression, na.rm = TRUE), .groups = "drop")

pdf("/mnt/nas-safu02/sdigiove_workspace/check_th_TCHRONET/new_tchroent_ties/ball/pictures/violin_trends.pdf"  ,
  height = 10, 
  width = 7,
  family = "ArialMT",
  useDingbats = FALSE)

ggplot(
  all_melted,
  aes(x = original_column, y = expression, fill = group)
) +
  geom_violin(
      trim = F,
      scale = "width",
      color = "black",
      alpha = 0.7
      #linewidth = 0.2
    ) +
  geom_boxplot(
      width = 0.15,
      outlier.shape = NA,
      color = "black",
      fill = "black",
      linewidth = 0.2
    ) +
  
  geom_hline(
    data = median_df,
    aes(yintercept = median_expr,),
    linetype = "dotted",
    linewidth = 0.8,
    color = "red"
  ) +
  facet_grid(community ~ group, scales = "free_x") +
  # facet_wrap(~community, scales = "free_x", ncol = 1) +
  scale_fill_manual(values = facet_colors) +
  # coord_cartesian(ylim = c(-5, 5)) +
  labs(y = "Z-score", x = NULL) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 15),
    axis.text.x = element_text(size = 15, angle = 90, vjust = 0.5),
    axis.text.y = element_text(size = 15),
    legend.position = "top",
    legend.text = element_text(size=15)

  )  +
  scale_x_discrete(
    labels = c(
      "patient_6_Healthy_ATACseq_REP1.mLb.clN.sorted.bam" = "H6",
 "patient_9_Healthy_ATACseq_REP1.mLb.clN.sorted.bam" = "H9",
 "patient_4_Healthy_ATACseq_REP1.mLb.clN.sorted.bam"    = "H4",
 "patient_7_Healthy_ATACseq_REP1.mLb.clN.sorted.bam" = "H7",
 "patient_5_Healthy_ATACseq_REP1.mLb.clN.sorted.bam" = "H5",
 "patient_8_Healthy_ATACseq_REP1.mLb.clN.sorted.bam" = "H8",
 "patient_14_Primary_ATACseq_REP1.mLb.clN.sorted.bam" = "P14",
 "patient_12_Primary_ATACseq_REP1.mLb.clN.sorted.bam" = "P12",
 "patient_10_Primary_ATACseq_REP1.mLb.clN.sorted.bam" = "P10",
 "patient_13_Primary_ATACseq_REP1.mLb.clN.sorted.bam" = "P13",
 "patient_11_Primary_ATACseq_REP1.mLb.clN.sorted.bam" = "P11",
 "patient_16_Primary_ATACseq_REP1.mLb.clN.sorted.bam" = "P16",
 "patient_7_Primary_ATACseq_REP1.mLb.clN.sorted.bam" = "P7",
 "patient_19_Primary_ATACseq_REP1.mLb.clN.sorted.bam" = "P19",
 "patient_17_Primary_ATACseq_REP1.mLb.clN.sorted.bam" = "P17",
 "patient_27_Primary_ATACseq_REP1.mLb.clN.sorted.bam" = "P27",
 "patient_26_Primary_ATACseq_REP1.mLb.clN.sorted.bam" = "P26",
 "patient_21_Remission_ATACseq_REP1.mLb.clN.sorted.bam" = "REM21",
 "patient_22_Remission_ATACseq_REP1.mLb.clN.sorted.bam" = "REM22",
 "patient_18_Remission_ATACseq_REP1.mLb.clN.sorted.bam" = "REM18",
 "patient_26_Remission_ATACseq_REP1.mLb.clN.sorted.bam" = "REM26",
 "patient_20_Remission_ATACseq_REP1.mLb.clN.sorted.bam" = "REM20",
 "patient_25_Remission_ATACseq_REP1.mLb.clN.sorted.bam" = "REM25",
 "patient_23_Remission_ATACseq_REP1.mLb.clN.sorted.bam" = "REM23",
 "patient_24_Remission_ATACseq_REP1.mLb.clN.sorted.bam" = "REM24",
 "patient_31_Relapse_ATACseq_REP1.mLb.clN.sorted.bam" = "REL31",
 "patient_28_Relapse_ATACseq_REP1.mLb.clN.sorted.bam" = "REL28",
 "patient_29_Relapse_ATACseq_REP1.mLb.clN.sorted.bam" = "REL29",
 "patient_33_Relapse_ATACseq_REP1.mLb.clN.sorted.bam" = "REL33",
 "patient_32_Relapse_ATACseq_REP1.mLb.clN.sorted.bam" = "REL32",
 "patient_27_Relapse_ATACseq_REP1.mLb.clN.sorted.bam" = "REL27",
 "patient_30_Relapse_ATACseq_REP1.mLb.clN.sorted.bam" = "REL30",
 "patient_14_Relapse_ATACseq_REP1.mLb.clN.sorted.bam"  = "REL14"
    ))
dev.off()



msigdbr_collections(db_species = "Hs")
genesets <- msigdbr(species = "Homo sapiens" , db_species = 'HS', collection = "C4" , subcollection = "3CA" ) #
genesets_removed  <- genesets |> dplyr::select(gs_name ,gene_symbol )
i = 1
j=0
for (x in ball_obj@GREAT_targets) {
    target_genes = unname(unlist(as.data.frame(x)$annotated_genes))
    x_enrichr <- enricher( target_genes , TERM2GENE = genesets_removed ) #   , universe = unname(unlist(reults[[1]]))
    x_df  <- x_enrichr@result |> filter(pvalue < 0.01) |> arrange( -FoldEnrichment )

    if(nrow(x_df) == 0) {    
        i = i + 1
        next
    }

    x_df  <- top_n(x_df,5,FoldEnrichment)
    x_df['community'] <- i

    if (j == 0) {
        final_df <- x_df
        j=1
    } else {
        final_df <- rbind(final_df , x_df)
    }
    i = i + 1
}

final_df_selected_Columns  <- final_df |> dplyr::select(community, Description , FoldEnrichment ,qvalue )
final_df_selected_Columns['qval_log10']  <- -log10(final_df_selected_Columns$qvalue)

clusteing_df = final_df_selected_Columns |> dplyr::select(Description , FoldEnrichment , community) |> spread(community , FoldEnrichment) |> column_to_rownames("Description")
clusteing_df[is.na(clusteing_df)] = 0
row.order <- stats::hclust(stats::dist(clusteing_df))$order
col.order <- stats::hclust(stats::dist(t(clusteing_df)))$order
clusteing_df[clusteing_df == 0] <-  NA
clusteing_df = clusteing_df[row.order, col.order]

final_df_selected_Columns  <- final_df_selected_Columns |> mutate(name = fct_relevel(Description, 
            rownames(clusteing_df)))
final_df_selected_Columns$community  <- factor(final_df_selected_Columns$community , levels = sort(final_df_selected_Columns |> pull(community) |> unique()))

final_df_selected_Columns$Description <- final_df_selected_Columns$Description |>
  str_remove("GAVISH_3CA_") |>
    str_replace_all("_", " ")|>
  str_replace_all("MALIGNANT", "MAL") |>
  str_replace_all("METAPROGRAM", "MP")


pdf("/mnt/nas-safu02/sdigiove_workspace/check_th_TCHRONET/new_tchroent_ties/ball/pictures/Metaprograms.pdf"  ,
  height = 5, 
  width = 9,
  family = "ArialMT",
  useDingbats = FALSE)
ggplot(final_df_selected_Columns , aes(x = community , y = Description , color = qval_log10 , size = FoldEnrichment)) +
  geom_point() + #color =Percent  , , size = 5
  scale_color_gradient(low = "blue" , high = "red" , na.value ="white")+
  scale_radius() +#trans = "log2"
  theme_classic() +
  # scale_y_discrete(limits=rownames(clusteing_df)) +
  theme(legend.position="top" ,  legend.box = "vertical" , legend.text = element_text(size=10)) +
  scale_size(range = c(3, 10)) +
  theme(text = element_text(size = 10))
dev.off()

### Save communities hg38
for (x in seq_along(list_communities_11)) {
  comm <- list_communities_11[[x]]
  comm_hg38 <- ball_obj@lifted_coords[ball_obj@lifted_coords$node %in% comm , ] |> pull(lifted_coord)
  comm_hg38 <- comm_hg38[!is.na(comm_hg38)]
  comm_hg38 <- comm_hg38 |> as.data.frame() |> separate(col = 1 , c("chr" , "start", "end") , sep = "-") |> dplyr::select(chr,start,end) 
  out_file <- paste("/mnt/nas-safu02/sdigiove_workspace/check_th_TCHRONET/new_tchroent_ties/ball/beds/community_hg38_" , as.character(x) , ".bed" , sep ="")
  write_delim(comm_hg38 , out_file , delim ="\t" , col_names = F)
}

pdf("/mnt/nas-safu02/sdigiove_workspace/check_th_TCHRONET/new_tchroent_ties/ball/pictures/TFs.pdf"  ,
  height = 8, 
  width = 7,
  family = "ArialMT",
  useDingbats = FALSE)
plot_cistrom("/mnt/nas-safu02/sdigiove_workspace/check_th_TCHRONET/new_tchroent_ties/ball/cistrome" , tf_name_file = "/mnt/nas-safu01/analysis/scripts/ScriptSdigiove/RegNetATACProject/T-ChroNet/paper_analysis/TFs_screening/TF_names_v_1.01.txt" )+
  theme(text = element_text(size = 15))+
  theme(legend.position="top" , legend.text = element_text(size=15))
dev.off()


pdf("/mnt/nas-safu02/sdigiove_workspace/check_th_TCHRONET/new_tchroent_ties/ball/pictures/Heatmap.pdf"  ,
  height = 9, 
  width = 5,
  family = "ArialMT",
  useDingbats = FALSE)
plot_heatmap_counts(ball_obj)
dev.off()



communities_11 <- ball_obj@clusters[,c('node' , "clusters_1")]
list_communities_11 <- list()
for (x in sort(unique(communities_11[["clusters_1"]])) ) {
  nodes <- communities_11[communities_11["clusters_1"] == x, ] |> pull(node)
  list_communities_11[[x]] <- nodes
}

communities_paper <- list.files("/mnt/nas-safu01/analysis/scripts/ScriptSdigiove/RegNetATACProject/T-ChroNet/paper_analysis/data/BCP-ALL/communities/" , full.names = T)
list_communities_paper <- list()
for (x in seq_along(communities_paper)) {
  #print(x)
  # print(communities_paper[[x]])
  list_communities_paper[[x]] <- read_delim(communities_paper[x], delim ="\t" , col_names = F) |> unite('nodes' , c(X1,X2,X3) , sep ="-") |> pull(nodes)
}


overlap_counts_sp_control <- overlap_ratios(
  list_communities_11,
  list_communities_paper
)

overlap_counts_sp_control <- t(overlap_counts_sp_control)
rownames(overlap_counts_sp_control) <- paste0('Community_old_' , c(1:4) , sep = "")
colnames(overlap_counts_sp_control) <- paste0('Community_' , c(1:4) , sep = "")


mat <- overlap_counts_sp_control

pdf("/mnt/nas-safu02/sdigiove_workspace/check_th_TCHRONET/new_tchroent_ties/ball/pictures/overlaps_communities_with_olds.pdf"  ,
  height = 6, 
  width = 6,
  family = "ArialMT",
  useDingbats = FALSE)
ht <- ComplexHeatmap::Heatmap(mat,
  name = "Overlaps",
cluster_rows = F,
cluster_columns = F,
show_column_dend = F,
show_row_dend = F ,
rect_gp = gpar(col = "white", lwd = 5),
cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(sprintf("%.2f", mat[i, j]), x, y, gp = gpar(fontsize = 10,    # change text size
        col = "white"))
  },
  col = viridis::viridis(100),
heatmap_legend_param = list(
    direction = "horizontal",
    title_position = "topcenter",
    legend_height = unit(8, "mm"),          # increases the color bar thickness
    title_gp = gpar(fontsize = 12, fontface = "bold"),  # legend title size
    labels_gp = gpar(fontsize = 10 , col="black")         
  ))
draw(ht, heatmap_legend_side = "top")
dev.off()
