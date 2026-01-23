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
  show_column_dend = F
  )
  
  return(h)
}

plot_cistrome <- function(
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


matrix_path <- "~/T-ChroNet/paper_analysis/data/LiverDevelopment/normalized_samplemean_multicov_all_sites_all_timepoint.tsv"
seiries_obj_liver <- build_TChroNetSeries_object("~/T-ChroNet/paper_analysis/data/LiverDevelopment/th",matrix_path,method = "Leiden",resolutions_list = seq(0.5, 1.5, 0.1),run_cd = T , transitivity = FALSE)


series_liver@best_th <- find_best_th(series_liver)

pdf("~/T-ChroNet/paper_analysis/data/LiverDevelopment/rand_index_map.pdf"  ,
  height = 9, 
  width = 5,
  family = "ArialMT",
  useDingbats = FALSE)
plot_randindex_map(series_liver )
dev.off()

pdf("~/T-ChroNet/paper_analysis/data/LiverDevelopment/sankey_plot_res_1.0.pdf"  ,
  height = 5, 
  width = 9,
  family = "ArialMT",
  useDingbats = FALSE)
plot_sankey_fixed_resolution(series_liver , resolution = 1.1)
dev.off()

pdf("~/T-ChroNet/paper_analysis/data/LiverDevelopment/relative_lcc.pdf"  ,
  height = 5, 
  width = 9,
  family = "ArialMT",
  useDingbats = FALSE)
plot_metrics(series_liver , "relative_lcc")
dev.off()


liver_obj <- build_TCrhoNetNetwork_from_series(series_liver , threshold = 0.8)

pdf("~/T-ChroNet/paper_analysis/data/LiverDevelopment/modularity.pdf"  ,
  height = 5, 
  width = 9,
  family = "ArialMT",
  useDingbats = FALSE)
plot_modularity(liver_obj)
dev.off()

pdf("~/T-ChroNet/paper_analysis/data/LiverDevelopment/community_sankey_th08.pdf"  ,
  height = 5, 
  width = 9,
  family = "ArialMT",
  useDingbats = FALSE)
plot_community_sankey(liver_obj , threshold = 0.8)
dev.off()

liver_obj <- find_best_resolution(liver_obj)

pdf("~/T-ChroNet/paper_analysis/data/LiverDevelopment/community_size.pdf"  ,
  height = 5, 
  width = 9,
  family = "ArialMT",
  useDingbats = FALSE)
liver_obj@clusters |> group_by(clusters_1.1) |> 
  summarise(peaks = n()) |> 
  mutate(clusters_1.1 = as.factor(clusters_1.1)) |> 
ggplot(aes(x = clusters_1.1 , y = peaks)) +
  geom_bar(stat = 'identity' , fill = "black") +
  theme_classic() + 
  coord_flip()+
  xlab("")+
  ylab("Size")+
  theme(text = element_text(size = 15))
dev.off()


liver_obj <- annotate_regions_from_bed(liver_obj , bed_path = "~/T-ChroNet/toy_data/data/mm10-cCREs.bed" , genome = 'mm10')

pdf("~/T-ChroNet/paper_analysis/data/LiverDevelopment/stachek_annotations_th08_res1.0.pdf"  ,
  height = 5, 
  width = 9,
  family = "ArialMT",
  useDingbats = FALSE)
plot_stacked_annotation(liver_obj , resolution = 1.1) +
  coord_flip() +
  theme(legend.position = "top" , text = element_text(size = 15)) +
  ylab("") +
  xlab("")
dev.off()




liver_obj <- run_GREAT_analysis(liver_obj, resolution = 1.1, genome = "mm10")


pdf("~/T-ChroNet/paper_analysis/data/LiverDevelopment/TFs.pdf"  ,
  height = 8, 
  width = 7,
  family = "ArialMT",
  useDingbats = FALSE)
plot_cistrome("~/T-ChroNet/paper_analysis/data/LiverDevelopment/cistrome/" , tf_name_file = "~/T-ChroNet/paper_analysis/TFs_screening/mouse_tfs.txt" )+
  theme(text = element_text(size = 15))
dev.off()

pdf("~/T-ChroNet/paper_analysis/data/LiverDevelopment/trends_violin.pdf"  ,
  height = 10, 
  width = 7,
  family = "ArialMT",
  useDingbats = FALSE)
plot_trends_zscore(liver_obj , resolution = 1.1)
dev.off()

pdf("~/T-ChroNet/paper_analysis/data/LiverDevelopment/trends_lines.pdf"  ,
  height = 10, 
  width = 7,
  family = "ArialMT",
  useDingbats = FALSE)
plot_line_trends_zscore(liver_obj , resolution = 1.1)
dev.off()

pdf("~/T-ChroNet/paper_analysis/data/LiverDevelopment/headtmaps_trends.pdf"  ,
  height = 9, 
  width = 5,
  family = "ArialMT",
  useDingbats = FALSE)
plot_heatmap_counts(liver_obj)
dev.off()

library(clusterProfiler)
library(msigdbr)
library(DOSE)

final_list_genes <- list()
for (x in liver_obj@GREAT_targets) {
  target_genes = unname(unlist(as.data.frame(x)$annotated_genes))
  final_list_genes <- c(final_list_genes , target_genes)
  final_list_genes <- unname(unlist(final_list_genes))
}

msigdbr_collections(db_species = "Mm") |> View()
genesets <- msigdbr(species = "Mus musculus" , db_species = 'MM', collection = "M8" )
genesets_removed  <- genesets |> dplyr::select(gs_name ,gene_symbol )
i = 1
j=0
for (x in liver_obj@GREAT_targets) {
  print(x)
    target_genes = unname(unlist(as.data.frame(x)$annotated_genes))
    x_enrichr <- enricher( target_genes , TERM2GENE = unique(genesets_removed) ) #   , universe = unname(unlist(reults[[1]])) #, universe = final_list_genes
    x_df  <- x_enrichr@result |> filter(qvalue < 0.05)
    x_df$FoldEnrichment <- (parse_ratio(x_df$GeneRatio)) /
                      (parse_ratio(x_df$BgRatio))
    x_df <- x_df |> arrange ( -FoldEnrichment )

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
  str_remove("DESCARTES_") |>
  str_remove("TABULA_MURIS") |>
  str_replace_all("_", " ")

pdf("~/T-ChroNet/paper_analysis/data/LiverDevelopment/Cell_annotations.pdf"  ,
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
  theme(legend.position="top" , legend.text = element_text(size=10))+
  theme(text = element_text(size=10))+
  xlab("")+
  ylab("")
  #scale_x_discrete(limits =colnames(clusteing_df), guide = guide_axis(angle = 90)) +
  # theme(axis.title.y=element_blank(), #axis.text.y = element_text(color = "grey20", size = 10, angle = 0, hjust = 1, vjust = 0, face = "plain") ,
  #       axis.text.y=element_blank(),
  #       legend.position="none")
dev.off()

### Save communities hg38
for (x in seq_along(list_communities_11)) {
  comm <- list_communities_11[[x]]
  comm_hg38 <- comm |> as.data.frame() |> separate(col = 1 , c("chr" , "start", "end") , sep = "-") |> dplyr::select(chr,start,end) 
  # comm_hg38 <- data.frame(list(peaks = comm))
  out_file <- paste("~/T-ChroNet/paper_analysis/data/LiverDevelopment/communities/community_mm10_" , as.character(x) , ".bed" , sep ="")
  write_delim(comm_hg38 , out_file , delim ="\t" , col_names = F)
}