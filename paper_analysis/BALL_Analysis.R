library(clustree)
library(tidyverse)
library(readr)
library(tidyverse)
library(ComplexHeatmap)
library(clusterProfiler)
library(msigdbr)


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
  show_column_dend = F
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
      select(Factor, max_value)

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

edge_files_ball <- list.files("/Users/sdigiove/Downloads/rds/th/" , full.names = T)
matrix_path <- "/Users/sdigiove/Downloads/rds/allPatientsPatientsNOMERGED_Ball_Multicov_nfcore.tsv"
series_ball_low <- build_TChroNetSeries_object(edge_files_ball,matrix_path,method = "Leiden",resolutions_list = seq(0.5, 1.5, 0.1),run_cd = T,seed = 123, transitivity = FALSE)
# write_rds(series_ball_low , "/mnt/nas-safu02/sdigiove_workspace/check_th_TCHRONET/ball/series_obj_ball_low.rds")

series_ball <- read_rds("/Users/sdigiove/Downloads/rds/series_obj_ball.rds")

ball_obj <- make_TChroNetNetwork_obj(edge_files = edge_files_ball , matrix_path = matrix_path , threshold = 0.9)

ball_obj <- compute_membership_nodes(ball_obj , resolutions = c(0.8 : 1.2))
# plot number of components and size distribution
m <- leidenAlg::find_partition_with_rep(ball_obj@graph, resolution = 0.8, edge_weights = E(ball_obj@graph)$weight , nrep = 5)

plot_randindex_map(series_ball)
plot_randindex_map_old(series_ball)


plot_metrics(series_ball , "relative_lcc") +
  geom_point( size = 3)+
  xlab("Threshold")+
  ylab("Relative LCC")+
  theme(text = element_text(size = 15))
ggsave("/Users/sdigiove/Documents/Work/Projects/T-ChroNet/picture_review/ball/relative_lcc.png" , height = 5 , width = 9 , units = 'in', dpi = 300)

series_ball@best_th <- find_best_th(series_ball) 

ball_obj <- build_TCrhoNetNetwork_from_series(series_ball )
write_rds(ball_obj , "/Users/sdigiove/Documents/Work/Projects/T-ChroNet/picture_review/ball/obj_ball.rds")

plot_community_sankey(ball_obj , threshold = 0.5)
ggsave("/Users/sdigiove/Documents/Work/Projects/T-ChroNet/picture_review/ball/community_sankey.png", height = 5 , width = 9 , units = 'in', dpi = 300)

plot_modularity(ball_obj)
ggsave("/Users/sdigiove/Documents/Work/Projects/T-ChroNet/picture_review/ball/modularity.png", height = 5 , width = 9 , units = 'in', dpi = 300)

ball_obj <- find_best_resolution(ball_obj)


plot_sankey_fixed_resolution(ball_obj)
ggsave("/Users/sdigiove/Documents/Work/Projects/T-ChroNet/picture_review/ball/sankey_plot_best_res.png" , height = 5 , width = 9 , units = 'in', dpi = 300)

ball_obj <- lift_network_coordinates(ball_obj , chain_path = "~/T-ChroNet/toy_data/data/hg19ToHg38.over.chain")
ball_obj <- annotate_regions_from_bed(liver_obj , bed_path = "~/T-ChroNet/toy_data/data/GRCh38-cCREs.bed" , genome = 'hg38')
plot_stacked_annotation(liver_obj , resolution = as.numeric(ball_obj@resolution)) +
  coord_flip() +
  theme(legend.position = "top" , text = element_text(size = 15)) +
  ylab("") +
  xlab("")
ggsave("/Users/sdigiove/Documents/Work/Projects/T-ChroNet/picture_review/ball/stachek_annotations_bestth_bestres.png", height = 5 , width = 9 , units = 'in', dpi = 300)

ball_obj <- run_GREAT_analysis(ball_obj, resolution = as.numeric(ball_obj@resolution), genome = "hg19")

plot_trends_zscore(ball_obj)+
  theme(text = element_text(size = 15))
ggsave("/Users/sdigiove/Documents/Work/Projects/T-ChroNet/picture_review/ball/boxplot_trends.png", height = 9 , width = 5 , units = 'in', dpi = 300)
plot_line_trends_zscore(ball_obj)+
    theme(text = element_text(size = 15))
ggsave("/Users/sdigiove/Documents/Work/Projects/T-ChroNet/picture_review/ball/lineplot_trends.png", height = 9 , width = 5 , units = 'in', dpi = 300)


msigdbr_collections(db_species = "Hs")
genesets <- msigdbr(species = "Homo sapiens" , db_species = 'HS', collection = "C4" , subcollection = "3CA" ) #
genesets_removed  <- genesets |> select(gs_name ,gene_symbol )
i = 1
j=0
for (x in ball_obj@GREAT_targets) {
    target_genes = unname(unlist(as.data.frame(x)$annotated_genes))
    x_enrichr <- enricher( target_genes , TERM2GENE = genesets_removed ) #   , universe = unname(unlist(reults[[1]]))
    x_df  <- x_enrichr@result |> filter(qvalue < 0.05) |> arrange( -FoldEnrichment )

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

final_df_selected_Columns  <- final_df |> select(community, Description , FoldEnrichment ,qvalue )
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

ggplot(final_df_selected_Columns , aes(x = community , y = Description , color = qval_log10 , size = FoldEnrichment)) +
  geom_point() + #color =Percent  , , size = 5
  scale_color_gradient(low = "blue" , high = "red" , na.value ="white")+
  scale_radius() +#trans = "log2"
  theme_classic() +
  # scale_y_discrete(limits=rownames(clusteing_df)) +
  theme(legend.position="top" , legend.text = element_text(size=10)) +
  scale_size(range = c(3, 10)) +
  theme(text = element_text(size = 15))
ggsave("/Users/sdigiove/Documents/Work/Projects/T-ChroNet/picture_review/ball/Metaprograms.png", height = 5 , width = 9 , units = 'in', dpi = 300)


best_clust_th <- paste("clusters_" , as.character(ball_obj@resolution) , sep ="")
ball_obj@clusters |> group_by(!!best_clust_th) |> 
  summarise(peaks = n()) |> 
  mutate(clusters_0.9 = as.factor(!!best_clust_th)) |> 
ggplot(aes(x = !!best_clust_th , y = peaks)) +
  geom_bar(stat = 'identity' , fill = "black") +
  theme_classic() + 
  coord_flip()+
  xlab("")+
  ylab("Size")+
  theme(text = element_text(size = 15))
ggsave("/Users/sdigiove/Documents/Work/Projects/T-ChroNet/picture_review/ball/community_size.png", height = 5 , width = 9 , units = 'in', dpi = 300)

