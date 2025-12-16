setwd("~/T-ChroNet/TChroNetR")
devtools::load_all()


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


library(clustree)
library(tidyverse)
library(readr)
library(tidyverse)
library(ComplexHeatmap)
library(clusterProfiler)
library(msigdbr)

# edge_files <- list.files("/mnt/nas-safu02/sdigiove_workspace/check_th_TCHRONET/thp1/th/" , full.names = T)
# matrix_path <- "/mnt/nas-safu01/analysis/scripts/ScriptSdigiove/RegNetATACProject/T-ChroNet/paper_analysis/data/THP1/lognorm_edgeR_limma_countsInCellReport_TheiCounts_NoStatic_mean.tsv"
# seiries_obj <- build_TChroNetSeries_object(edge_files,matrix_path ,method = "Leiden",resolutions_list = seq(0.5, 1.5, 0.1),run_cd = T,seed = 1234 , transitivity = FALSE)
# write_rds(seiries_obj , "/mnt/nas-safu02/sdigiove_workspace/check_th_TCHRONET/thp1/series_obj_thp1.rds")

# edge_files_liver <- list.files("/mnt/nas-safu02/sdigiove_workspace/check_th_TCHRONET/liver/th/" , full.names = T)
# matrix_path <- "/mnt/nas-safu01/analysis/scripts/ScriptSdigiove/RegNetATACProject/T-ChroNet/paper_analysis/data/LiverDevelopment/normalized_samplemean_multicov_all_sites_all_timepoint.tsv"
# seiries_obj_liver <- build_TChroNetSeries_object(edge_files_liver,matrix_path,method = "Leiden",resolutions_list = seq(0.5, 1.5, 0.1),run_cd = T,seed = 123, transitivity = FALSE)
# write_rds(seiries_obj_liver , "/mnt/nas-safu02/sdigiove_workspace/check_th_TCHRONET/liver/series_obj_liver.rds")

# edge_files_ball <- list.files("/mnt/nas-safu02/sdigiove_workspace/check_th_TCHRONET/ball/th/" , full.names = T)
# matrix_path <- "/mnt/nas-safu01/analysis/scripts/ScriptSdigiove/RegNetATACProject/T-ChroNet/paper_analysis/data/BCP-ALL/allPatientsPatientsNOMERGED_Ball_Multicov_nfcore.tsv"
# series_ball <- build_TChroNetSeries_object(edge_files_ball,matrix_path,method = "Leiden",resolutions_list = seq(0.5, 1.5, 0.1),run_cd = T,seed = 123, transitivity = FALSE)
# write_rds(series_ball , "/mnt/nas-safu02/sdigiove_workspace/check_th_TCHRONET/ball/series_obj_ball.rds")


#### THP1
series_thp1 <- read_rds("/Users/sdigiove/Downloads/rds/series_obj_thp1.rds")

plot_randindex_map(series_thp1 ,"/mnt/nas-safu02/sdigiove_workspace/check_th_TCHRONET/thp1/pictures/rand_index_map.png" )
plot_sankey_fixed_resolution(series_thp1 , resolution = 0.9)
ggsave("/Users/sdigiove/Downloads/rds/thp1/pictures/sankey_plot_res_1_5.png" , height = 5 , width = 9 , units = 'in', dpi = 300)
plot_metrics(series_thp1 , "relative_lcc") +
  theme(text = element_text(size = 15))
#ggsave("/mnt/nas-safu02/sdigiove_workspace/check_th_TCHRONET/thp1/pictures/relative_lcc.png" , height = 5 , width = 9 , units = 'in', dpi = 300)

series_thp1 <- find_best_th(series_thp1)

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

### Liver
series_liver <- read_rds("/Users/sdigiove/Downloads/rds/series_obj_liver.rds")

plot_randindex_map(series_liver ,"/Users/sdigiove/Documents/Work/Projects/T-ChroNet/picture_review/liver/rand_index_map.png" )
plot_sankey_fixed_resolution(series_liver , resolution = 0.8)
ggsave("/Users/sdigiove/Documents/Work/Projects/T-ChroNet/picture_review/liver/sankey_plot_res_1.1.png" , height = 5 , width = 9 , units = 'in', dpi = 300)
plot_metrics(series_liver , "relative_lcc")
ggsave("/Users/sdigiove/Documents/Work/Projects/T-ChroNet/picture_review/liver/relative_lcc.png" , height = 5 , width = 9 , units = 'in', dpi = 300)

liver_obj <- build_TCrhoNetNetwork_from_series(series_liver , threshold = 0.8)
plot_modularity(liver_obj)
ggsave("/Users/sdigiove/Documents/Work/Projects/T-ChroNet/picture_review/liver/modularity.png", height = 5 , width = 9 , units = 'in', dpi = 300)
plot_community_sankey(liver_obj , threshold = 0.8)
ggsave("/Users/sdigiove/Documents/Work/Projects/T-ChroNet/picture_review/liver/community_sankey_th08.png", height = 5 , width = 9 , units = 'in', dpi = 300)

liver_obj <- find_best_resolution(liver_obj)

liver_obj@clusters |> group_by(clusters_1.1) |> 
  summarise(peaks = n()) |> 
  mutate(clusters_0.9 = as.factor(clusters_1.1)) |> 
ggplot(aes(x = clusters_1.1 , y = peaks)) +
  geom_bar(stat = 'identity' , fill = "black") +
  theme_classic() + 
  coord_flip()+
  xlab("")+
  ylab("Size")+
  theme(text = element_text(size = 15))
ggsave("/Users/sdigiove/Documents/Work/Projects/T-ChroNet/picture_review/liver/community_size.png", height = 5 , width = 9 , units = 'in', dpi = 300)


liver_obj <- annotate_regions_from_bed(liver_obj , bed_path = "~/T-ChroNet/toy_data/data/mm10-cCREs_mod.bed" , genome = 'mm10')
plot_stacked_annotation(liver_obj , resolution = 1.1) +
  coord_flip() +
  theme(legend.position = "top" , text = element_text(size = 15)) +
  ylab("") +
  xlab("")
ggsave("/Users/sdigiove/Documents/Work/Projects/T-ChroNet/picture_review/liver/stachek_annotations_th08_res1.1.png", height = 5 , width = 9 , units = 'in', dpi = 300)


communities_08 <- liver_obj@clusters[,c('node' , "clusters_0.8")]
list_communities_08 <- list()
for (x in sort(unique(communities_08[["clusters_0.8"]])) ) {
  nodes <- communities_08[communities_08["clusters_0.8"] == x, ] |> pull(node)
  list_communities_08[[x]] <- nodes
}

communities_11 <- liver_obj@clusters[,c('node' , "clusters_1.1")]
list_communities_11 <- list()
for (x in sort(unique(communities_11[["clusters_1.1"]])) ) {
  nodes <- communities_11[communities_11["clusters_1.1"] == x, ] |> pull(node)
  list_communities_11[[x]] <- nodes
}

communities_paper <- list.files("~/T-ChroNet/paper_analysis/data/LiverDevelopment/communities/" , full.names = T)
list_communities_paper <- list()
for (x in seq_along(communities_paper)) {
  #print(x)
  # print(communities_paper[[x]])
  list_communities_paper[[x]] <- read_delim(communities_paper[x], delim ="\t" , col_names = F) |> unite('nodes' , c(X1,X2,X3) , sep ="-") |> pull(nodes)
}

listA <- list_communities_11
names(listA) <- paste0("TChroNetR_" , as.character(c(1:5)))
listB <- list_communities_08
names(listB) <- paste0("paper_" , as.character(c(1:4)))

intersections <- expand.grid(A = names(listA), B = names(listB))

intersections$overlap <- mapply(
  function(a, b) length(intersect(listA[[a]], listB[[b]])),
  intersections$A,
  intersections$B
)
library(UpSetR)

# Get universe of elements
all_elems <- unique(unlist(c(listA, listB)))

# Convert each list into binary membership matrix
make_binary <- function(lst) {
  sapply(lst, function(x) as.numeric(all_elems %in% x))
}

binA <- make_binary(listA)
binB <- make_binary(listB)

# Combine sets
binary_matrix <- cbind(binA, binB)

# png("/mnt/nas-safu02/sdigiove_workspace/check_th_TCHRONET/thp1/pictures/paper_vs_TCNR.png", height = 5 , width = 9 , units = 'in', res = 300)
upset(
  as.data.frame(binary_matrix),
  sets = colnames(binary_matrix),
  keep.order = T
  # order.by = "freq",
)
# dev.off()

liver_obj <- run_GREAT_analysis(liver_obj, resolution = 1.1, genome = "mm10")

### Save communities hg38
for (x in seq_along(list_communities_11)) {
  comm <- list_communities_11[[x]]
  comm_hg38 <- comm |> as.data.frame() |> separate(col = 1 , c("chr" , "start", "end") , sep = "-") |> select(chr,start,end) 
  out_file <- paste("/Users/sdigiove/Documents/Work/Projects/T-ChroNet/picture_review/liver/bed/community_hg38_" , as.character(x) , ".bed" , sep ="")
  write_delim(comm_hg38 , out_file , delim ="\t" , col_names = F)
}

plot_cistrom("/Users/sdigiove/Downloads/cistrom_liver/" , tf_name_file = "~/T-ChroNet/paper_analysis/TFs_screening/mouse_tfs.txt" )+
  theme(text = element_text(size = 15))
ggsave("/Users/sdigiove/Documents/Work/Projects/T-ChroNet/picture_review/liver/TFs.png", width = 7, height = 8 , units = 'in', dpi = 300)


plot_trends_zscore(liver_obj , resolution = 1.1)
# ggsave("/Users/sdigiove/Documents/Work/Projects/T-ChroNet/picture_review/liver/trends_violin.png", height = 9 , width = 5 , units = 'in', dpi = 300)
plot_line_trends_zscore(liver_obj , resolution = 1.1)
# ggsave("/Users/sdigiove/Documents/Work/Projects/T-ChroNet/picture_review/liver/trends_lines.png", height = 9 , width = 5 , units = 'in', dpi = 300)

# png("/Users/sdigiove/Documents/Work/Projects/T-ChroNet/picture_review/liver/Heatmap_trends.png", height = 9 , width = 5 , units = 'in', res = 300)
plot_heatmap_counts(liver_obj)
# dev.off()

library(clusterProfiler)
library(msigdbr)

final_list_genes <- list()
for (x in liver_obj@GREAT_targets) {
  target_genes = unname(unlist(as.data.frame(x)$annotated_genes))
  final_list_genes <- c(final_list_genes , target_genes)
  final_list_genes <- unname(unlist(final_list_genes))
}

msigdbr_collections(db_species = "Mm") |> View()
genesets <- msigdbr(species = "Mus musculus" , db_species = 'MM', collection = "M8" )
genesets_removed  <- genesets |> select(gs_name ,gene_symbol )
i = 1
j=0
for (x in liver_obj@GREAT_targets) {
    target_genes = unname(unlist(as.data.frame(x)$annotated_genes))
    x_enrichr <- enricher( target_genes , TERM2GENE = genesets_removed  , universe = final_list_genes) #   , universe = unname(unlist(reults[[1]]))
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
  str_remove("DESCARTES_") |>
  str_remove("TABULA_MURIS") |>
  str_replace_all("_", " ")

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
ggsave("/Users/sdigiove/Documents/Work/Projects/T-ChroNet/picture_review/liver/Cell_annotations.png", height = 5 , width = 9 , units = 'in', dpi = 300)


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
  


### B-ALL
series_ball <- read_rds("/Users/sdigiove/Downloads/rds/series_obj_ball.rds")

plot_randindex_map(series_ball)

plot_metrics(series_ball , "relative_lcc")
ggsave("/Users/sdigiove/Documents/Work/Projects/T-ChroNet/picture_review/ball/relative_lcc.png" , height = 5 , width = 9 , units = 'in', dpi = 300)

series_ball@best_th <- find_best_th(series_ball) 


series_ball <- build_TCrhoNetNetwork_from_series(series_ball )

plot_community_sankey(series_ball , threshold = 0.8)
ggsave("/Users/sdigiove/Documents/Work/Projects/T-ChroNet/picture_review/ball/community_sankey_th08.png", height = 5 , width = 9 , units = 'in', dpi = 300)

plot_modularity(series_ball)
ggsave("/Users/sdigiove/Documents/Work/Projects/T-ChroNet/picture_review/ball/modularity.png", height = 5 , width = 9 , units = 'in', dpi = 300)
series_ball <- find_best_resolution(series_ball)
plot_sankey_fixed_resolution(series_ball , resolution = 0.8)
ggsave("/Users/sdigiove/Documents/Work/Projects/T-ChroNet/picture_review/ball/sankey_plot_res_1.1.png" , height = 5 , width = 9 , units = 'in', dpi = 300)

liver_obj@clusters |> group_by(clusters_1.1) |> 
  summarise(peaks = n()) |> 
  mutate(clusters_0.9 = as.factor(clusters_1.1)) |> 
ggplot(aes(x = clusters_1.1 , y = peaks)) +
  geom_bar(stat = 'identity' , fill = "black") +
  theme_classic() + 
  coord_flip()+
  xlab("")+
  ylab("Size")+
  theme(text = element_text(size = 15))
ggsave("/Users/sdigiove/Documents/Work/Projects/T-ChroNet/picture_review/liver/community_size.png", height = 5 , width = 9 , units = 'in', dpi = 300)


liver_obj <- annotate_regions_from_bed(liver_obj , bed_path = "~/T-ChroNet/toy_data/data/mm10-cCREs_mod.bed" , genome = 'mm10')
plot_stacked_annotation(liver_obj , resolution = 1.1) +
  coord_flip() +
  theme(legend.position = "top" , text = element_text(size = 15)) +
  ylab("") +
  xlab("")
ggsave("/Users/sdigiove/Documents/Work/Projects/T-ChroNet/picture_review/liver/stachek_annotations_th08_res1.1.png", height = 5 , width = 9 , units = 'in', dpi = 300)


communities_08 <- liver_obj@clusters[,c('node' , "clusters_0.8")]
list_communities_08 <- list()
for (x in sort(unique(communities_08[["clusters_0.8"]])) ) {
  nodes <- communities_08[communities_08["clusters_0.8"] == x, ] |> pull(node)
  list_communities_08[[x]] <- nodes
}

communities_11 <- liver_obj@clusters[,c('node' , "clusters_1.1")]
list_communities_11 <- list()
for (x in sort(unique(communities_11[["clusters_1.1"]])) ) {
  nodes <- communities_11[communities_11["clusters_1.1"] == x, ] |> pull(node)
  list_communities_11[[x]] <- nodes
}

communities_paper <- list.files("~/T-ChroNet/paper_analysis/data/LiverDevelopment/communities/" , full.names = T)
list_communities_paper <- list()
for (x in seq_along(communities_paper)) {
  #print(x)
  # print(communities_paper[[x]])
  list_communities_paper[[x]] <- read_delim(communities_paper[x], delim ="\t" , col_names = F) |> unite('nodes' , c(X1,X2,X3) , sep ="-") |> pull(nodes)
}

listA <- list_communities_11
names(listA) <- paste0("TChroNetR_" , as.character(c(1:5)))
listB <- list_communities_08
names(listB) <- paste0("paper_" , as.character(c(1:4)))

intersections <- expand.grid(A = names(listA), B = names(listB))

intersections$overlap <- mapply(
  function(a, b) length(intersect(listA[[a]], listB[[b]])),
  intersections$A,
  intersections$B
)
library(UpSetR)

# Get universe of elements
all_elems <- unique(unlist(c(listA, listB)))

# Convert each list into binary membership matrix
make_binary <- function(lst) {
  sapply(lst, function(x) as.numeric(all_elems %in% x))
}

binA <- make_binary(listA)
binB <- make_binary(listB)

# Combine sets
binary_matrix <- cbind(binA, binB)

# png("/mnt/nas-safu02/sdigiove_workspace/check_th_TCHRONET/thp1/pictures/paper_vs_TCNR.png", height = 5 , width = 9 , units = 'in', res = 300)
upset(
  as.data.frame(binary_matrix),
  sets = colnames(binary_matrix),
  keep.order = T
  # order.by = "freq",
)
# dev.off()

liver_obj <- run_GREAT_analysis(liver_obj, resolution = 1.1, genome = "mm10")

### Save communities hg38
for (x in seq_along(list_communities_11)) {
  comm <- list_communities_11[[x]]
  comm_hg38 <- comm |> as.data.frame() |> separate(col = 1 , c("chr" , "start", "end") , sep = "-") |> select(chr,start,end) 
  out_file <- paste("/Users/sdigiove/Documents/Work/Projects/T-ChroNet/picture_review/liver/bed/community_hg38_" , as.character(x) , ".bed" , sep ="")
  write_delim(comm_hg38 , out_file , delim ="\t" , col_names = F)
}

library(tidyverse)
library(ComplexHeatmap)

# THP-1
folders_to_compare = list.files('/Users/sdigiove/Downloads/cistrom_liver/')

for (community_num in seq_along(folders_to_compare) ){ 
    
    community <- folders_to_compare[[community_num]]
    to_read = paste('/Users/sdigiove/Downloads/cistrom_liver/' ,community , sep = "")

    homer_tfs = read_delim(to_read , delim = ",")

    #homer_tfs  <- homer_tfs |> filter(grepl('liver' , Biosource ) | grepl('Liver' , Biosource ))
    homer_tfs = homer_tfs |> group_by(Factor) |> summarise(max_value = max(GIGGLE_score)) |> select(Factor , max_value)
    

    colnames(homer_tfs) = c('Factor' , community_num)

    if (community == folders_to_compare[1]){
        final_homer_tfs = homer_tfs}
    else {
        final_homer_tfs = merge.data.frame(final_homer_tfs , homer_tfs , all = T)}
}
final_homer_tfs <- final_homer_tfs |> column_to_rownames("Factor")

tfs_name = read_delim("~/T-ChroNet/paper_analysis/TFs_screening/TF_names_v_1.01.txt" , delim = "\t" , col_names = F)
tfs_name <- tfs_name %>%
  mutate(across(where(is.character), toupper))
real_tfs  <- c()
for (x in rownames(final_homer_tfs)){
    if (x %in% tfs_name[['X1']]){
        real_tfs <- c(real_tfs , x)
    }
}
final_homer_tfs_only_tfs <- final_homer_tfs[rownames(final_homer_tfs) %in% real_tfs,]
saved_tfs = c()
for (column_name in colnames(final_homer_tfs_only_tfs)) {
    final_homer_tfs_dropped  <- final_homer_tfs_only_tfs[!is.na(final_homer_tfs_only_tfs[column_name]), column_name, drop = FALSE]  
    final_homer_tfs_dropped <- final_homer_tfs_dropped[order(-final_homer_tfs_dropped[[column_name]]), column_name, drop = FALSE]
    top_10 = head(final_homer_tfs_dropped, 10)
    saved_tfs <- c(saved_tfs, rownames(top_10))
}
saved_tfs <- unique(saved_tfs)

cluster_map <- final_homer_tfs_only_tfs[saved_tfs,]
cluster_map[is.na(cluster_map)] <- 0
r.cluster = hclust(stats::dist(cluster_map, method = "euclidean"), method="ward.D2")

final_homer_scatterplot <- final_homer_tfs_only_tfs[saved_tfs,]|> rownames_to_column("Factor") |> gather(value = 'ComboScore' , key = 'Community' , -Factor) 
final_homer_scatterplot$Factor <- factor(final_homer_scatterplot$Factor, levels = r.cluster$labels[r.cluster$order])
final_homer_scatterplot$Community <- as.numeric(gsub("\\D", "", final_homer_scatterplot$Community))

ggplot(data = final_homer_scatterplot , aes(x = Community, y = Factor ) ) +
geom_point( aes( size = ComboScore) ,  shape = 21  , stroke = 0.6  , color = 'black') +
geom_point( aes( size = ComboScore) ,  shape = 21  , fill = "black" , stroke = 0.6 , alpha= 0.4 , color = 'black') +
theme_classic() +
ylab("")+
theme(legend.position="top" , legend.text = element_text(size=10) )+
  theme(text = element_text(size = 15))
  #axis.text.y = element_text(color = "grey20", size = 10, angle = 0, hjust = 1, vjust = 0) ) +
#scale_x_discrete( limits = unique(final_homer_scatterplot$Community) )
# theme(axis.title.y=element_blank(), #axis.text.y = element_text(color = "grey20", size = 10, angle = 0, hjust = 1, vjust = 0, face = "plain") ,
#       axis.text.y=element_blank(),
#       legend.position="none")
ggsave("/Users/sdigiove/Documents/Work/Projects/T-ChroNet/picture_review/liver/TFs.png", width = 7, height = 8 , units = 'in', dpi = 300)


plot_trends_zscore(liver_obj , resolution = 1.1)
# ggsave("/Users/sdigiove/Documents/Work/Projects/T-ChroNet/picture_review/liver/trends_violin.png", height = 9 , width = 5 , units = 'in', dpi = 300)
plot_line_trends_zscore(liver_obj , resolution = 1.1)
# ggsave("/Users/sdigiove/Documents/Work/Projects/T-ChroNet/picture_review/liver/trends_lines.png", height = 9 , width = 5 , units = 'in', dpi = 300)

# png("/Users/sdigiove/Documents/Work/Projects/T-ChroNet/picture_review/liver/Heatmap_trends.png", height = 9 , width = 5 , units = 'in', res = 300)
plot_heatmap_counts(liver_obj)
# dev.off()

library(clusterProfiler)
library(msigdbr)

final_list_genes <- list()
for (x in liver_obj@GREAT_targets) {
  target_genes = unname(unlist(as.data.frame(x)$annotated_genes))
  final_list_genes <- c(final_list_genes , target_genes)
  final_list_genes <- unname(unlist(final_list_genes))
}

msigdbr_collections(db_species = "Mm") |> View()
genesets <- msigdbr(species = "Mus musculus" , db_species = 'MM', collection = "M8" )
genesets_removed  <- genesets |> select(gs_name ,gene_symbol )
i = 1
j=0
for (x in liver_obj@GREAT_targets) {
    target_genes = unname(unlist(as.data.frame(x)$annotated_genes))
    x_enrichr <- enricher( target_genes , TERM2GENE = genesets_removed  , universe = final_list_genes) #   , universe = unname(unlist(reults[[1]]))
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
  str_remove("DESCARTES_") |>
  str_remove("TABULA_MURIS") |>
  str_replace_all("_", " ")

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
ggsave("/Users/sdigiove/Documents/Work/Projects/T-ChroNet/picture_review/liver/Cell_annotations.png", height = 5 , width = 9 , units = 'in', dpi = 300)


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

object <- series_ball
color_limits = c(0.0, 1)
