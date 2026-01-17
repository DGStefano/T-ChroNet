library(tidyverse)
library(readr)
library(tidyverse)
library(ComplexHeatmap)
library(clusterProfiler)
library(msigdbr)
library(TChroNetR)
library(dplyr)

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


##### LOADING TCHRONET FILES AND CREATE THE SERIES OBJECT

# edge_files <- list.files("/mnt/nas-safu02/sdigiove_workspace/check_th_TCHRONET/new_tchroent_ties/thp1/th/" , full.names = T)
# matrix_path <- "~/T-ChroNet/paper_analysis/data/THP1/lognorm_edgeR_limma_countsInCellReport_TheiCounts_NoStatic_mean.tsv"
# series_thp1 <- build_TChroNetSeries_object(edge_files,matrix_path ,method = "Leiden",resolutions_list = seq(0.5, 1.5, 0.1),run_cd = T,seed = 1234 , transitivity = FALSE)
# write_rds(series_thp1 , "/mnt/nas-safu02/sdigiove_workspace/check_th_TCHRONET/new_tchroent_ties/thp1/series_obj_thp1.rds")

##### READING THE SERIES OBJECT RDS
series_thp1 <- read_rds("/mnt/nas-safu02/sdigiove_workspace/check_th_TCHRONET/new_tchroent_ties/thp1/series_obj_thp1.rds")

##### FINDING THE BEST THRESHOLD
series_thp1@best_th <- find_best_th(series_thp1)

##### PLOTTING THE RAND INDEX MATRIX OF NEIGHBOURS THRESHOLD AND RESOLUTIONS
pdf("/mnt/nas-safu02/sdigiove_workspace/check_th_TCHRONET/new_tchroent_ties/thp1/pictures/pdf/rand_index_map.pdf" ,
  height = 10, 
  width = 8,
  family = "ArialMT",
  useDingbats = FALSE)
plot_randindex_map(series_thp1)
dev.off()

##### PLOTTING THE COMMUNITIES STABILITY ALONG THRESHOLDS
pdf("/mnt/nas-safu02/sdigiove_workspace/check_th_TCHRONET/new_tchroent_ties/thp1/pictures/pdf/sankey_plot_res_1_3.pdf" ,
  height = 5, 
  width = 9,
  family = "ArialMT",
  useDingbats = FALSE)
plot_sankey_fixed_resolution(series_thp1 , resolution = 1.0)
dev.off()

##### PLOTTING THE REALTIVE LARGEST CONNECTED COMPONENTS FOR EACH THRESHOLD
to_plot <- plot_metrics(series_thp1 , "relative_lcc") +
  theme(text = element_text(size = 15))

pdf("/mnt/nas-safu02/sdigiove_workspace/check_th_TCHRONET/new_tchroent_ties/thp1/pictures/pdf/relative_lcc.pdf" ,
  height = 5, 
  width = 9,
  family = "ArialMT",
  useDingbats = FALSE)
print(to_plot)
dev.off()

##### BUILDING TCHRONET OBJECT AT RESOLUTION OF 0.8 STARTING FROM THE SERIES
thp1_th08 <- build_TCrhoNetNetwork_from_series(series_thp1 , threshold = 0.8)
# write_rds(thp1_th08 , "/mnt/nas-safu02/sdigiove_workspace/check_th_TCHRONET/new_tchroent_ties/thp1/Network_obj_thp1.rds")
# thp1_th08 <- read_rds("/mnt/nas-safu02/sdigiove_workspace/check_th_TCHRONET/new_tchroent_ties/thp1/Network_obj_thp1.rds")

##### PLOT MODULARITY FOR ALL THE RESOLUTIONS
to_plot <- plot_modularity(thp1_th08)+
  theme(text = element_text(size = 15))
pdf("/mnt/nas-safu02/sdigiove_workspace/check_th_TCHRONET/new_tchroent_ties/thp1/pictures/pdf/modularity_th08.pdf" ,
  height = 5, 
  width = 9,
  family = "ArialMT",
  useDingbats = FALSE)
print(to_plot)
dev.off()

to_plot <- plot_community_sankey(thp1_th08 , threshold = 0.8)
pdf("/mnt/nas-safu02/sdigiove_workspace/check_th_TCHRONET/new_tchroent_ties/thp1/pictures/pdf/community_sankey_th08.pdf" ,
  height = 5, 
  width = 9,
  family = "ArialMT",
  useDingbats = FALSE)
print(to_plot)
dev.off()

thp1_th08 <- find_best_resolution(thp1_th08)

to_plot <- thp1_th08@clusters |> group_by(clusters_1) |> 
  summarise(peaks = n()) |> 
  mutate(clusters_1 = as.factor(clusters_1)) |> 
ggplot(aes(x = clusters_1 , y = peaks)) +
  geom_bar(stat = 'identity' , fill = "black") +
  theme_classic() + 
  coord_flip()+
  xlab("")+
  ylab("Size")+
  theme(text = element_text(size = 15))
pdf("/mnt/nas-safu02/sdigiove_workspace/check_th_TCHRONET/new_tchroent_ties/thp1/pictures/pdf/community_size.pdf" ,
  height = 5, 
  width = 9,
  family = "ArialMT",
  useDingbats = FALSE)
print(to_plot)
dev.off()

thp1_th08 <- lift_network_coordinates(thp1_th08 , chain_path = "~/T-ChroNet/toy_data/data/hg19ToHg38.over.chain")
thp1_th08 <- annotate_regions_from_bed(thp1_th08 , bed_path = "~/T-ChroNet/toy_data/data/GRCh38-cCREs.bed" , genome = "hg38")
to_plot <- plot_stacked_annotation(thp1_th08) +
  coord_flip() +
  theme(legend.position = "top" , text = element_text(size = 15)) +
  ylab("") +
  xlab("")
pdf("/mnt/nas-safu02/sdigiove_workspace/check_th_TCHRONET/new_tchroent_ties/thp1/pictures/pdf/stachek_annotations_th08_res1_0.pdf" ,
  height = 5, 
  width = 9,
  family = "ArialMT",
  useDingbats = FALSE)
print(to_plot)
dev.off()


thp1_th08 <- run_GREAT_analysis(thp1_th08, genome = "hg19")

communities_09 <- thp1_th08@clusters[,c('node' , "clusters_1")]
list_communities_09 <- list()
for (x in sort(unique(communities_09[["clusters_1"]])) ) {
  nodes <- communities_09[communities_09["clusters_1"] == x, ] |> pull(node)
  list_communities_09[[x]] <- nodes
}

### Save communities hg38
for (x in seq_along(list_communities_09)) {
  comm <- list_communities_09[[x]]
  comm_hg38 <- thp1_th08@lifted_coords |> dplyr::filter(node %in% comm ) |> dplyr::select(chr_hg38,start_hg38,end_hg38) 
  out_file <- paste("/mnt/nas-safu02/sdigiove_workspace/check_th_TCHRONET/new_tchroent_ties/thp1/beds/hg38/community_hg38_" , as.character(x) , ".bed" , sep ="")
  # comm_hg38 <- data.frame(list(peaks = comm))
  #out_file <- paste("/mnt/nas-safu02/sdigiove_workspace/check_th_TCHRONET/new_tchroent_ties/thp1/beds/community_hg19_" , as.character(x) , ".bed" , sep ="")
  write_delim(comm_hg38 , out_file , delim ="\t" , col_names = F)
}


to_plot <- plot_trends_zscore(thp1_th08)+
  theme(text = element_text(size = 15))

pdf("/mnt/nas-safu02/sdigiove_workspace/check_th_TCHRONET/new_tchroent_ties/thp1/pictures/pdf/boxplot_trends.pdf" ,
  height = 9, 
  width = 5,
  family = "ArialMT",
  useDingbats = FALSE)
print(to_plot)
dev.off()


to_plot <- plot_line_trends_zscore(thp1_th08)+
    theme(text = element_text(size = 15))

pdf("/mnt/nas-safu02/sdigiove_workspace/check_th_TCHRONET/new_tchroent_ties/thp1/pictures/pdf/lineplot_trends.pdf" ,
  height = 9, 
  width = 5,
  family = "ArialMT",
  useDingbats = FALSE)
print(to_plot)
dev.off()


msigdbr_collections(db_species = "Hs")
genesets <- msigdbr(species = "Homo sapiens" , db_species = 'HS', collection = "H") #, subcollection = "CP:REACTOME"
genesets_removed  <- genesets |> dplyr::select(gs_name ,gene_symbol )
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
  str_remove("HALLMARK_") |>
  str_replace_all("_", " ")

to_plot <- ggplot(final_df_selected_Columns , aes(x = community , y = Description , color = qval_log10 , size = FoldEnrichment)) +
  geom_point() + #color =Percent  , , size = 5
  scale_color_gradient(low = "blue" , high = "red" , na.value ="white")+
  scale_radius() +#trans = "log2"
  theme_classic() +
  # scale_y_discrete(limits=rownames(clusteing_df)) +
  theme(legend.position="top" , legend.text = element_text(size=10)) +
  scale_size(range = c(3, 10)) +
  theme(text = element_text(size = 15))
pdf("/mnt/nas-safu02/sdigiove_workspace/check_th_TCHRONET/new_tchroent_ties/thp1/pictures/pdf/Hallmarks.pdf",
  height = 5, 
  width = 9,
  family = "ArialMT",
  useDingbats = FALSE)
print(to_plot)
dev.off()



to_plot <- plot_cistrom("/mnt/nas-safu02/sdigiove_workspace/check_th_TCHRONET/new_tchroent_ties/thp1/beds/cistrome/" , tf_name_file = "~/T-ChroNet/paper_analysis/TFs_screening/TF_names_v_1.01.txt", ) +
  theme(text = element_text(size = 15))
pdf( "/mnt/nas-safu02/sdigiove_workspace/check_th_TCHRONET/new_tchroent_ties/thp1/pictures/pdf/TFs.pdf",
  height = 8, 
  width = 7,
  family = "ArialMT",
  useDingbats = FALSE)
print(to_plot)
dev.off()




pdf("/mnt/nas-safu02/sdigiove_workspace/check_th_TCHRONET/new_tchroent_ties/thp1/pictures/pdf/heatmap.pdf",
  height = 9, 
  width = 5,
  family = "ArialMT",
  useDingbats = FALSE)
plot_heatmap_counts(thp1_th08)
dev.off()


communities_11 <- thp1_th08@clusters[,c('node' , "clusters_1")]
list_communities_11 <- list()
for (x in sort(unique(communities_11[["clusters_1"]])) ) {
  nodes <- communities_11[communities_11["clusters_1"] == x, ] |> pull(node)
  list_communities_11[[x]] <- nodes
}

genes_df <- read_delim("~/T-ChroNet/paper_analysis/data/THP1/Genes_lognorm_edgeR_limma_countsInCellReport_TheiCounts_All_mean.tsv", delim ="\t" ) |> column_to_rownames('hgnc')

all_long <- map_dfr(seq_along(list_communities_11), function(i) {

  # ------------------------------
  # Community name
  # ------------------------------
  comm_name <- paste0("Community_", i)

  # ------------------------------
  # ATAC peaks
  # ------------------------------
  peaks <- list_communities_11[[i]] #unname(unlist(thp1_th08@GREAT_targets[[i]]@))[unname(unlist(thp1_th08@GREAT_targets[[i]]$dist_to_TSS)) < 10000 & unname(unlist(thp1_th08@GREAT_targets[[i]]$dist_to_TSS)) > -10000] # list_communities_11[[i]]
  peaks_use <- intersect(peaks, rownames(thp1_th08@matrix))

  atac_mat <- thp1_th08@matrix[peaks_use, , drop = FALSE]

  atac_z <- t(scale(t(atac_mat)))

  atac_long <- as.data.frame(atac_z) |>
    rownames_to_column("feature") |>
    pivot_longer(
      cols = -feature,
      names_to = "time",
      values_to = "value"
    ) |>
    mutate(
      Dataset = "ATAC",
      community = comm_name
    )

  # ------------------------------
  # RNA genes
  # ------------------------------
  genes <- unname(unlist(thp1_th08@GREAT_targets[[i]]$annotated_genes))
  genes_use <- intersect(genes, rownames(genes_df))

  rna_mat <- genes_df[genes_use, , drop = FALSE]

  rna_z <- t(scale(t(rna_mat)))

  rna_long <- as.data.frame(rna_z) |>
    rownames_to_column("feature") |>
    pivot_longer(
      cols = -feature,
      names_to = "time",
      values_to = "value"
    ) |>
    mutate(
      Dataset = "RNA",
      community = comm_name
    )

  # ------------------------------
  # Combine
  # ------------------------------
  bind_rows(atac_long, rna_long)
})
time_levels <- colnames(thp1_th08@matrix)

time_labels <- c(
  "0min","30min","60min","90min",
  "120min","240min","360min","1440min"
)

all_long$time <- factor(all_long$time, levels = time_levels)

cols <- c("ATAC" = "orange", "RNA" = "blue")

p <- ggplot(
  all_long,
  aes(x = time, y = value, color = Dataset, group = Dataset)
) +
  scale_color_manual(values = cols)+
  stat_summary(
    fun = mean,
    geom = "line",
    linewidth = 1
  ) +
  stat_summary(
    fun = mean,
    geom = "point",
    size = 2
  ) +
  facet_wrap(
     ~ community,
    scales = "free_y"
  ) +
  scale_x_discrete(labels = time_labels) +
  labs(
    x = NULL,
    y = "Z-score"
  ) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 9),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top",
    text = element_text(size = 15)
  )
pdf("/mnt/nas-safu02/sdigiove_workspace/check_th_TCHRONET/new_tchroent_ties/thp1/pictures/pdf/rna_genes.pdf" ,
  height = 8, 
  width = 7,
  family = "ArialMT",
  useDingbats = FALSE)
print(p)
dev.off()



sites_th1_cellreport <- read_delim(
  "~/T-ChroNet/paper_analysis/data/THP1/their_regions_timepoints.bed", delim = "\t"
) 

sites_th1_cellreport <- sites_th1_cellreport |>
  mutate(
    sites = paste(chr, start, end, sep = "-")
  ) |>
  dplyr::select(sites, cluster) |>
  filter(sites %in% rownames(thp1_th08@matrix))
sites_th1_cellreport[!duplicated(sites_th1_cellreport$sites),]

communities_name <- c(
  "dn.early","dn.grad","dn.mid","dn.var",
  "up.grad","up.late","up.mid","up.var"
)

communities_list_comparison <- lapply(
  communities_name,
  function(x) {
    sites_th1_cellreport |>
      filter(cluster == x) |>
      pull(sites)
  }
)

count_overlaps <- function(list1, list2) {
  # list1, list2: lists of vectors (character or numeric)

  overlaps <- vector("list", length(list1))

  for (i in seq_along(list1)) {
    overlap_row <- numeric(length(list2))

    for (j in seq_along(list2)) {
      common_elements <- intersect(list1[[i]], list2[[j]])
      overlap_row[j] <- length(common_elements)
    }

    overlaps[[i]] <- overlap_row
  }

  # Convert to matrix for easier downstream use
  do.call(rbind, overlaps)
}
overlap_ratios <- function(list1, list2) {
  # list1, list2: lists of vectors
  # Returns matrix with rows = list1, cols = list2

  ratios <- vector("list", length(list1))

  for (i in seq_along(list1)) {
    ratio_row <- numeric(length(list2))

    for (j in seq_along(list2)) {
      common_elements <- intersect(list1[[i]], list2[[j]])

      if (length(list2[[j]]) > 0) {
        ratio_row[j] <- length(common_elements) / length(list2[[j]])
      } else {
        ratio_row[j] <- 0
      }
    }

    ratios[[i]] <- ratio_row
  }

  do.call(rbind, ratios)
}



overlap_counts_sp_control <- overlap_ratios(
  communities_list_comparison,
  list_communities_11
)

overlap_counts_sp_control <- t(overlap_counts_sp_control)
colnames(overlap_counts_sp_control) <- communities_name
rownames(overlap_counts_sp_control) <- paste0('Community_' , c(1:5) , sep = "")


mat <- overlap_counts_sp_control

pdf("/mnt/nas-safu02/sdigiove_workspace/check_th_TCHRONET/new_tchroent_ties/thp1/pictures/overlaps_communities.pdf"  ,
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
    labels_gp = gpar(fontsize = 10)         
  ))
draw(ht, heatmap_legend_side = "top")
dev.off()


#### Overlaps old new

old_comms <- list()
files_old_list <- list.files(path = "/mnt/nas-safu02/sdigiove_workspace/check_th_TCHRONET/thp1/beds/" , full.names = TRUE)
for(x in seq_along(files_old_list)) {
  df <- read_delim(files_old_list[[x]] , delim = "\t" , col_names = FALSE)
  old_comms[[x]] <- df |> pull(X1)
}


overlap_counts_sp_control <- overlap_ratios(
  list_communities_11,
  old_comms
)

overlap_counts_sp_control <- t(overlap_counts_sp_control)
rownames(overlap_counts_sp_control) <- paste0('Community_old_' , c(1:6) , sep = "")
colnames(overlap_counts_sp_control) <- paste0('Community_' , c(1:5) , sep = "")


mat <- overlap_counts_sp_control

pdf("/mnt/nas-safu02/sdigiove_workspace/check_th_TCHRONET/new_tchroent_ties/thp1/pictures/overlaps_communities_with_olds.pdf"  ,
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


##### Making table of communities and original annotation
tchronet_communities <- thp1_th08@clusters |> dplyr::select(node , clusters_1)

merged_annotations_clusters <- merge.data.frame(sites_th1_cellreport , tchronet_communities , by.x = 'sites', by.y = 'node' , all = T)

merged_annotations_clusters |> separate(sites , c('chr','start','end') , sep ="-") |> write_delim( "/mnt/nas-safu02/sdigiove_workspace/check_th_TCHRONET/new_tchroent_ties/thp1/annotations_and_clusters_thp1.tsv" , delim = "\t" )

matrix_path <- "~/T-ChroNet/paper_analysis/data/THP1/lognorm_edgeR_limma_countsInCellReport_TheiCounts_NoStatic_mean.tsv"
matrix_values <- read_delim(matrix_path , delim = "\t")
merged_annotations_clusters <- merge.data.frame(merged_annotations_clusters , matrix_values , by.x = 'sites',by.y = 'peaks')

merged_annotations_clusters |> separate(sites , c('chr','start','end') , sep ="-") |> write_delim( "/mnt/nas-safu02/sdigiove_workspace/check_th_TCHRONET/new_tchroent_ties/thp1/annotations_and_clusters_thp1.tsv" , delim = "\t" )
