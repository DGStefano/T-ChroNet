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

# edge_files <- list.files("/mnt/nas-safu02/sdigiove_workspace/check_th_TCHRONET/thp1/th/" , full.names = T)
# matrix_path <- "/mnt/nas-safu01/analysis/scripts/ScriptSdigiove/RegNetATACProject/T-ChroNet/paper_analysis/data/THP1/lognorm_edgeR_limma_countsInCellReport_TheiCounts_NoStatic_mean.tsv"
# seiries_obj <- build_TChroNetSeries_object(edge_files,matrix_path ,method = "Leiden",resolutions_list = seq(0.5, 1.5, 0.1),run_cd = T,seed = 1234 , transitivity = FALSE)
# write_rds(seiries_obj , "/mnt/nas-safu02/sdigiove_workspace/check_th_TCHRONET/thp1/series_obj_thp1.rds")



#### THP1
series_thp1 <- read_rds("/Users/sdigiove/Downloads/rds/series_obj_thp1.rds")

series_thp1 <- find_best_th(series_thp1)

plot_randindex_map(series_thp1 ,"/mnt/nas-safu02/sdigiove_workspace/check_th_TCHRONET/thp1/pictures/rand_index_map.png" )
plot_sankey_fixed_resolution(series_thp1 , resolution = 0.9)
ggsave("/Users/sdigiove/Downloads/rds/thp1/pictures/sankey_plot_res_1_5.png" , height = 5 , width = 9 , units = 'in', dpi = 300)
plot_metrics(series_thp1 , "relative_lcc") +
  theme(text = element_text(size = 15))
#ggsave("/mnt/nas-safu02/sdigiove_workspace/check_th_TCHRONET/thp1/pictures/relative_lcc.png" , height = 5 , width = 9 , units = 'in', dpi = 300)


thp1_th08 <- build_TCrhoNetNetwork_from_series(series_thp1 , threshold = 0.8)
plot_modularity(thp1_th08)+
  theme(text = element_text(size = 15))
# ggsave("/mnt/nas-safu02/sdigiove_workspace/check_th_TCHRONET/thp1/pictures/modularity_th08.png" , height = 5 , width = 9 , units = 'in', dpi = 300)
plot_community_sankey(thp1_th08 , threshold = 0.8)
# ggsave("/mnt/nas-safu02/sdigiove_workspace/check_th_TCHRONET/thp1/pictures/community_sankey_th08.png", height = 5 , width = 9 , units = 'in', dpi = 300)

thp1_th08 <- find_best_resolution(thp1_th08)
thp1_th08@clusters |> group_by(clusters_0.9) |> 
  summarise(peaks = n()) |> 
  mutate(clusters_0.9 = as.factor(clusters_0.9)) |> 
ggplot(aes(x = clusters_0.9 , y = peaks)) +
  geom_bar(stat = 'identity' , fill = "black") +
  theme_classic() + 
  coord_flip()+
  xlab("")+
  ylab("Size")+
  theme(text = element_text(size = 15))
# ggsave("/mnt/nas-safu02/sdigiove_workspace/check_th_TCHRONET/thp1/pictures/community_size.png", height = 5 , width = 9 , units = 'in', dpi = 300)

plot_sankey_fixed_resolution(series_thp1 , resolution = 0.9)
# ggsave("/mnt/nas-safu02/sdigiove_workspace/check_th_TCHRONET/thp1/pictures/sankey_fixed_resolution.png", height = 5 , width = 9 , units = 'in', dpi = 300)


thp1_th08 <- lift_network_coordinates(thp1_th08 , chain_path = "../toy_data/data/hg19ToHg38.over.chain")
thp1_th08 <- annotate_regions_from_bed(thp1_th08 , bed_path = "../toy_data/data/GRCh38-cCREs_2.bed" , genome = "hg38")
plot_stacked_annotation(thp1_th08) +
  coord_flip() +
  theme(legend.position = "top" , text = element_text(size = 15)) +
  ylab("") +
  xlab("")
# ggsave("/mnt/nas-safu02/sdigiove_workspace/check_th_TCHRONET/thp1/pictures/stachek_annotations_th08_res09.png", height = 5 , width = 9 , units = 'in', dpi = 300)

thp1_th08 <- run_GREAT_analysis(thp1_th08, genome = "hg19")

communities_09 <- thp1_th08@clusters[,c('node' , "clusters_0.9")]
list_communities_09 <- list()
for (x in sort(unique(communities_09[["clusters_0.9"]])) ) {
  nodes <- communities_09[communities_09["clusters_0.9"] == x, ] |> pull(node)
  list_communities_09[[x]] <- nodes
}

### Save communities hg38
for (x in seq_along(list_communities_09)) {
  comm <- list_communities_09[[x]]
  comm_hg38 <- thp1_th08@lifted_coords |> filter(node %in% comm ) |> select(chr_hg38,start_hg38,end_hg38) 
  out_file <- paste("/mnt/nas-safu02/sdigiove_workspace/check_th_TCHRONET/review_results/thp1/bed/community_hg38_" , as.character(x) , ".bed" , sep ="")
  write_delim(comm_hg38 , out_file , delim ="\t" , col_names = F)
}


plot_trends_zscore(thp1_th08)+
  theme(text = element_text(size = 15))
# ggsave("/mnt/nas-safu02/sdigiove_workspace/check_th_TCHRONET/thp1/pictures/boxplot_trends.png", height = 9 , width = 5 , units = 'in', dpi = 300)
plot_line_trends_zscore(thp1_th08)+
    theme(text = element_text(size = 15))
# ggsave("/mnt/nas-safu02/sdigiove_workspace/check_th_TCHRONET/thp1/pictures/lineplot_trends.png", height = 9 , width = 5 , units = 'in', dpi = 300)



msigdbr_collections(db_species = "Hs")
genesets <- msigdbr(species = "Homo sapiens" , db_species = 'HS', collection = "H") #, subcollection = "CP:REACTOME"
genesets_removed  <- genesets |> select(gs_name ,gene_symbol )
i = 1
j=0
for (x in thp1_th08@GREAT_targets) {
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
  str_remove("HALLMARK_") |>
  str_replace_all("_", " ")

ggplot(final_df_selected_Columns , aes(x = community , y = Description , color = qval_log10 , size = FoldEnrichment)) +
  geom_point() + #color =Percent  , , size = 5
  scale_color_gradient(low = "blue" , high = "red" , na.value ="white")+
  scale_radius() +#trans = "log2"
  theme_classic() +
  # scale_y_discrete(limits=rownames(clusteing_df)) +
  theme(legend.position="top" , legend.text = element_text(size=10)) +
  scale_size(range = c(3, 10)) +
  theme(text = element_text(size = 15))
ggsave("/mnt/nas-safu02/sdigiove_workspace/check_th_TCHRONET/thp1/pictures/Hallmarks.png", height = 5 , width = 9 , units = 'in', dpi = 300)

plot_cistrom("/Users/sdigiove/Downloads/cistromdb/" , tf_name_file = "~/T-ChroNet/paper_analysis/TFs_screening/TF_names_v_1.01.txt", ) +
  theme(text = element_text(size = 15))
ggsave("/Users/sdigiove/Documents/Work/Projects/T-ChroNet/picture_review/thp1/TFs_Cistrome.png",width = 7, height = 8 , units = 'in', dpi = 300)


png("/Users/sdigiove/Documents/Work/Projects/T-ChroNet/picture_review/thp1/heatmap.png", height = 9 , width = 5 , units = 'in', res = 300)
plot_heatmap_counts(thp1_th08)
dev.off()