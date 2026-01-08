library(WGCNA)
library(tidyverse)
library(flashClust)

allowWGCNAThreads(nThreads = 40)
powers = c(c(1:10), seq(from = 12, to = 40, by = 2))

input_mat <- read_tsv("/mnt/nas-safu01/analysis/scripts/ScriptSdigiove/RegNetATACProject/T-ChroNet/paper_analysis/data/LiverDevelopment/normalized_samplemean_multicov_all_sites_all_timepoint.tsv") |> tibble::column_to_rownames('peaks')
input_mat <- t(input_mat)

sft <- pickSoftThreshold(
  input_mat,
  powerVector = powers,
  verbose = 5,
  moreNetworkConcepts = TRUE
)

sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;

plot(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R^2",
     main = paste("Scale independence")
)
text(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = cex1, col = "red"
)
abline(h = 0.90, col = "red")
plot(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     type = "n",
     main = paste("Mean connectivity")
)
text(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     labels = powers,
     cex = cex1, col = "red")



softPower <- 1
adjacency <- WGCNA::adjacency(input_mat, power = softPower)
TOM=TOMsimilarityFromExpr(input_mat,networkType = "unsigned", TOMType = "unsigned", power = softPower);
dissTOM=1-TOM
geneTree = flashClust(as.dist(dissTOM),method="ward");

set.seed(123)
dynamicMods_100 = cutreeDynamic(dendro = geneTree,  method="tree", minClusterSize = 100);
set.seed(123)
dynamicMods_1000 = cutreeDynamic(dendro = geneTree,  method="tree", minClusterSize = 1000);

dynamicColors = labels2colors(dynamicMods)

plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "Gene dendrogram and module colors")

peak_clusters_100 <- data.frame(list(peaks = colnames(input_mat) , clusuters = dynamicMods_100))
colnames(peak_clusters_100) <- c('peaks' , "clusters_100")
write_delim(peak_clusters_100 , "/mnt/nas-safu02/sdigiove_workspace/check_th_TCHRONET/WGCNA/clusters_min100.tsv" , delim ="\t")
peak_clusters_1000 <- data.frame(list(peaks = colnames(input_mat) , clusuters = dynamicMods_1000))
colnames(peak_clusters_1000) <- c('peaks' , "clusters_1000")
write_delim(peak_clusters_1000 , "/mnt/nas-safu02/sdigiove_workspace/check_th_TCHRONET/WGCNA/clusters_min1000.tsv" , delim ="\t")


files_lists <- list.files("/mnt/nas-safu02/sdigiove_workspace/check_th_TCHRONET/liver/bed/" , full.names = T)
i = 0
for ( x in  seq_along(files_lists)) {
  df <- read_delim(files_lists[[x]] , delim = "\t" , col_names = F)
  df['cluster_TCNR'] <- x
  if ( i == 0 ){
    final_tcnr_clst <- df
    i = 1 
  }
  else {
    final_tcnr_clst <- rbind(final_tcnr_clst , df)
  }
}

final_comparison <- merge.data.frame(peak_clusters_100 , final_tcnr_clst , by.x = 'peaks' , by.y = 'X1')
final_comparison <- merge.data.frame(final_comparison , peak_clusters_1000 , by = 'peaks')


final_comparison$clusters_100 <- factor(final_comparison$clusters_100 , levels = sort(unique(final_comparison$clusters_100)))
final_comparison$cluster_TCNR <- factor(final_comparison$cluster_TCNR , levels = sort(unique(final_comparison$cluster_TCNR)))
final_comparison$clusters_1000 <- factor(final_comparison$clusters_1000 , levels = sort(unique(final_comparison$clusters_1000)))

library(ggplot2)
library(ggalluvial)

# Plotting the transition between 'clusters' and 'clusters_tcnr'
ggplot(final_comparison, aes(axis1 = clusters_100, axis2 = cluster_TCNR , axis3 = clusters_1000)) +
  geom_alluvium(aes(fill = cluster_TCNR), width = 1/12) +
  geom_stratum(width = 1/12, fill = "grey80", color = "black") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("WGCNA Clusters min 100", "TCNR Clusters" ,"WGCNA Clusters min 1000"), expand = c(.05, .05 , .05)) +
  labs(title = "Peaks Cluster Membership Transition",
       y = "Number of Peaks") +
  theme_minimal() +
  theme(text = element_text(size = 15))
ggsave("/mnt/nas-safu02/sdigiove_workspace/check_th_TCHRONET/WGCNA/picture/WGCNA_seq_100_our_1000_LIVER.pdf" , height = 5 , width = 9 , units = 'in' , dpi = 300)

library(rGREAT)
library(GenomicRanges)
interrogate_GREAT_regions <- function(gr,
                                      comm_name,
                                      species) {
  suppressPackageStartupMessages({
    library(rGREAT)
    library(dplyr)
  })
  
  job <- submitGreatJob(
    gr,
    version = "4.0.4",
    species = species,
    rule = "basalPlusExt",
    adv_upstream = 2.0,
    adv_downstream = 1.0,
    request_interval = 12
  )
  
  tbl <- getEnrichmentTables(job)
  
  target_genes <- getRegionGeneAssociations(job)
  
  return(target_genes)
}

final_target_genes <- list()  
for (comm in sort(unique(peak_clusters_100$clusters_100))) {
  
  regions <- peak_clusters_100[peak_clusters_100$clusters_100 == comm, ] |> select(peaks) |> separate(peaks,c('chr','start','end') , sep ="-")

  if (nrow(regions) == 0) next
    
  gr <- makeGRangesFromDataFrame(
    regions,
    seqnames.field = colnames(regions)[1],
    start.field = colnames(regions)[2],
    end.field = colnames(regions)[3],
    keep.extra.columns = TRUE
  )
    
  GREAT_results <- interrogate_GREAT_regions(
    gr = gr,
    comm_name = paste0("community_", comm),
    species = 'mm10'
  )
    
  
  final_target_genes[[comm +1 ]] <- unlist(unname(GREAT_results$annotated_genes))
}



library(clusterProfiler)
library(msigdbr)
library(DOSE)

msigdbr_collections(db_species = "Mm") |> View()
genesets <- msigdbr(species = "Mus musculus" , db_species = 'MM', collection = "M8" )
genesets_removed  <- genesets |> select(gs_name ,gene_symbol )
i = 1
j=0
for (x in final_target_genes) {
    target_genes = x
    x_enrichr <- enricher( target_genes , TERM2GENE = genesets_removed  ) #   , universe = unname(unlist(reults[[1]]))
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
ggsave("/mnt/nas-safu02/sdigiove_workspace/check_th_TCHRONET/WGCNA/picture/Cell_annotations_WGCNA_100_2.pdf",device = cairo_pdf , 
  height = 5, 
  width = 9, 
  units = 'in',
  dpi = 300
)

#### LIVER
powers = c(c(1:10), seq(from = 12, to = 40, by = 2))

input_mat <- read_tsv("/mnt/nas-safu01/analysis/scripts/ScriptSdigiove/RegNetATACProject/T-ChroNet/paper_analysis/data/LiverDevelopment/normalized_samplemean_multicov_all_sites_all_timepoint.tsv") |> tibble::column_to_rownames('peaks')input_mat <- read_tsv("/mnt/nas-safu01/analysis/scripts/ScriptSdigiove/RegNetATACProject/T-ChroNet/paper_analysis/data/THP1/lognorm_edgeR_limma_countsInCellReport_TheiCounts_NoStatic_mean.tsv") |> tibble::column_to_rownames('peaks')
input_mat <- t(input_mat)

sft <- pickSoftThreshold(
  input_mat,
  powerVector = powers,
  verbose = 5,
  moreNetworkConcepts = TRUE
)

sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;

plot(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R^2",
     main = paste("Scale independence")
)
text(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = cex1, col = "red"
)
abline(h = 0.90, col = "red")
plot(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     type = "n",
     main = paste("Mean connectivity")
)
text(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     labels = powers,
     cex = cex1, col = "red")



softPower <- 1
adjacency <- WGCNA::adjacency(input_mat, power = softPower)
TOM=TOMsimilarityFromExpr(input_mat,networkType = "unsigned", TOMType = "unsigned", power = softPower);
dissTOM=1-TOM
geneTree = flashClust(as.dist(dissTOM),method="ward");

dynamicMods_100 = cutreeDynamic(dendro = geneTree,  method="tree", minClusterSize = 100);
dynamicMods_1000 = cutreeDynamic(dendro = geneTree,  method="tree", minClusterSize = 1000);

dynamicColors = labels2colors(dynamicMods)

plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "Gene dendrogram and module colors")

peak_clusters_100 <- data.frame(list(peaks = colnames(input_mat) , clusuters = dynamicMods_100))
colnames(peak_clusters_100) <- c('peaks' , "clusters_100")
write_delim(peak_clusters_100 , "/mnt/nas-safu02/sdigiove_workspace/check_th_TCHRONET/WGCNA/clusters_min100_THP1.tsv" , delim ="\t")
peak_clusters_100 <- read_delim("/mnt/nas-safu02/sdigiove_workspace/check_th_TCHRONET/WGCNA/clusters_min100.tsv" , delim ="\t")
peak_clusters_1000 <- data.frame(list(peaks = colnames(input_mat) , clusuters = dynamicMods_1000))
colnames(peak_clusters_1000) <- c('peaks' , "clusters_1000")
write_delim(peak_clusters_1000 , "/mnt/nas-safu02/sdigiove_workspace/check_th_TCHRONET/WGCNA/clusters_min1000_THP1.tsv" , delim ="\t")
peak_clusters_1000 <- read_delim("/mnt/nas-safu02/sdigiove_workspace/check_th_TCHRONET/WGCNA/clusters_min1000.tsv" , delim ="\t")

files_lists <- list.files("/mnt/nas-safu02/sdigiove_workspace/check_th_TCHRONET/thp1/beds/" , full.names = T)
i = 0
for ( x in  seq_along(files_lists)) {
  df <- read_delim(files_lists[[x]] , delim = "\t" , col_names = F)
  df['cluster_TCNR'] <- x
  if ( i == 0 ){
    final_tcnr_clst <- df
    i = 1 
  }
  else {
    final_tcnr_clst <- rbind(final_tcnr_clst , df)
  }
}

final_comparison <- merge.data.frame(peak_clusters_100 , final_tcnr_clst , by.x = 'peaks' , by.y = 'X1')
final_comparison <- merge.data.frame(final_comparison , peak_clusters_1000 , by = 'peaks')


final_comparison$clusters_100 <- factor(final_comparison$clusters_100 , levels = sort(unique(final_comparison$clusters_100)))
final_comparison$cluster_TCNR <- factor(final_comparison$cluster_TCNR , levels = sort(unique(final_comparison$cluster_TCNR)))
final_comparison$clusters_1000 <- factor(final_comparison$clusters_1000 , levels = sort(unique(final_comparison$clusters_1000)))

library(ggplot2)
library(ggalluvial)

# Plotting the transition between 'clusters' and 'clusters_tcnr'
ggplot(final_comparison, aes(axis1 = clusters_100, axis2 = cluster_TCNR , axis3 = clusters_1000)) +
  geom_alluvium(aes(fill = cluster_TCNR), width = 1/12) +
  geom_stratum(width = 1/12, fill = "grey80", color = "black") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("WGCNA Clusters min 100", "TCNR Clusters" ,"WGCNA Clusters min 1000"), expand = c(.05, .05 , .05)) +
  labs(title = "Peaks Cluster Membership Transition",
       y = "Number of Peaks") +
  theme_minimal() +
  theme(text = element_text(size = 15))
ggsave("/mnt/nas-safu02/sdigiove_workspace/check_th_TCHRONET/WGCNA/picture/WGCNA_seq_100_our_1000_THP1.pdf" , height = 5 , width = 9 , units = 'in' , dpi = 300)


#### MonaLisa evaluation
library(monaLisa)
library(JASPAR2024)
library(TFBSTools)
JASPAR2024 <- JASPAR2024()
JASPARConnect <- RSQLite::dbConnect(RSQLite::SQLite(), db(JASPAR2024))
pwms <- getMatrixSet(JASPARConnect,
                     opts = list(
                      tax_group = "vertebrates",
                      collection="CORE",
                      matrixtype = "PWM"))

# disconnect Db
RSQLite::dbDisconnect(JASPARConnect)

genome_obj <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
monaLisa_enrichments <- list()
  for (i in sort(unique(peak_clusters_1000$clusters_1000))) {
    community_nodes <- peak_clusters_1000[peak_clusters_1000$clusters_1000 == i, ] |> dplyr::select(peaks) |> tidyr::separate(peaks,c('chr','start','end') , sep ="-")

    gr <- GenomicRanges::makeGRangesFromDataFrame(
    community_nodes,
    seqnames.field = colnames(community_nodes)[1],
    start.field = colnames(community_nodes)[2],
    end.field = colnames(community_nodes)[3],
    keep.extra.columns = TRUE
    )

    seqs <- BSgenome::getSeq(genome_obj, gr)

    se_genome <- monaLisa::calcBinnedMotifEnrR(
    seqs = seqs,
    pwmL = pwms,                      
    background = "genome",            
    genome = genome_obj, 
    genome.oversample = 2,
    BPPARAM = BiocParallel::MulticoreParam(20)
    )
    i = i+1
    monaLisa_enrichments[[i]] <- se_genome
}

negLog10Padj_th  = 4.0
log2enr_th = 1.0
top_n = 20

for (enrichments_num in seq_along(monaLisa_enrichments)) {
    enrichments <- monaLisa_enrichments[[enrichments_num]]
      # select strongly enriched motifs
    sel <- apply(SummarizedExperiment::assay(enrichments, "negLog10Padj"), 1,
                function(x) max(abs(x), 0, na.rm = TRUE)) > negLog10Padj_th
    sel_enr <- apply(SummarizedExperiment::assay(enrichments, "log2enr"), 1, 
                  function(x) max(x, na.rm = TRUE)) > log2enr_th # 2-fold enriched
    final_sel <- sel & sel_enr
    seSel <- enrichments[final_sel,]


    matrix_enrich <- seSel@assays@data$log2enr
    id_to_names <- seSel@elementMetadata[c("motif.id" , "motif.name")]
    rownames(matrix_enrich) <- unlist(unname(id_to_names[id_to_names$motif.id %in% rownames(matrix_enrich), "motif.name" ]))
    colnames(matrix_enrich) <- paste('community_' , as.character(enrichments_num) , sep ="")
    matrix_enrich <- matrix_enrich |> as.data.frame() |> rownames_to_column('tfs')
      
    if (enrichments_num == 1 ) {
      final_enrich <- matrix_enrich
    }
    else {
      final_enrich <- merge.data.frame(final_enrich , matrix_enrich , by = "tfs" , all = T)
    }
    }
    final_enrich_2 <- final_enrich |> column_to_rownames('tfs')


    saved_tfs <- c()

    for (column_name in colnames(final_enrich_2)) {
      col_df <- final_enrich_2[!is.na(final_enrich_2[[column_name]]), column_name, drop = FALSE]
      col_df <- col_df[order(-col_df[[column_name]]),, drop = FALSE]
      top_values <- head(col_df, top_n)
      saved_tfs <- c(saved_tfs, rownames(top_values))
    }

    saved_tfs <- unique(saved_tfs)

    # matrix for ordering
    cluster_map <- final_enrich_2[saved_tfs, ]
    cluster_map[is.na(cluster_map)] <- 0

    r.cluster <- fastcluster::hclust(dist(cluster_map, method = "euclidean"), method = "ward.D2")

    # -------------------------
    # 5. Long-format for plotting
    # -------------------------
    final_homer_scatterplot <- final_enrich_2[saved_tfs, ] |>
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
ggplot(final_homer_scatterplot, aes(x = Community, y = Factor)) +
      geom_point(aes(size = ComboScore),
                shape = 21, stroke = 0.6, color = "black") +
      geom_point(aes(size = ComboScore),
                shape = 21, fill = "black", stroke = 0.6, alpha = 0.4) +
      scale_x_discrete(drop=F)+
      theme_classic() +
      ylab("")+
  theme(text = element_text(size=10)) +
    theme(
      legend.position = "top",
      legend.text = element_text(size = 10)
    )
ggsave("/mnt/nas-safu02/sdigiove_workspace/check_th_TCHRONET/WGCNA/picture/monalisa_min1000.pdf",device = cairo_pdf , 
  height = 9, 
  width = 5, 
  units = 'in',
  dpi = 300
)



monaLisa_enrichments <- list()
  for (i in sort(unique(peak_clusters_100$clusters_100))) {
    community_nodes <- peak_clusters_100[peak_clusters_100$clusters_100 == i, ] |> dplyr::select(peaks) |> tidyr::separate(peaks,c('chr','start','end') , sep ="-")

    gr <- GenomicRanges::makeGRangesFromDataFrame(
    community_nodes,
    seqnames.field = colnames(community_nodes)[1],
    start.field = colnames(community_nodes)[2],
    end.field = colnames(community_nodes)[3],
    keep.extra.columns = TRUE
    )

    seqs <- BSgenome::getSeq(genome_obj, gr)

    se_genome <- monaLisa::calcBinnedMotifEnrR(
    seqs = seqs,
    pwmL = pwms,                      
    background = "genome",            
    genome = genome_obj, 
    genome.oversample = 2,
    BPPARAM = BiocParallel::MulticoreParam(20)
    )
    i = i +1
    monaLisa_enrichments[[i]] <- se_genome
}

negLog10Padj_th  = 4.0
log2enr_th = 1.0
top_n = 20

for (enrichments_num in seq_along(monaLisa_enrichments)) {
    enrichments <- monaLisa_enrichments[[enrichments_num]]
      # select strongly enriched motifs
    sel <- apply(SummarizedExperiment::assay(enrichments, "negLog10Padj"), 1,
                function(x) max(abs(x), 0, na.rm = TRUE)) > negLog10Padj_th
    sel_enr <- apply(SummarizedExperiment::assay(enrichments, "log2enr"), 1, 
                  function(x) max(x, na.rm = TRUE)) > log2enr_th # 2-fold enriched
    final_sel <- sel & sel_enr
    seSel <- enrichments[final_sel,]


    matrix_enrich <- seSel@assays@data$log2enr
    id_to_names <- seSel@elementMetadata[c("motif.id" , "motif.name")]
    rownames(matrix_enrich) <- unlist(unname(id_to_names[id_to_names$motif.id %in% rownames(matrix_enrich), "motif.name" ]))
    colnames(matrix_enrich) <- paste('community_' , as.character(enrichments_num) , sep ="")
    matrix_enrich <- matrix_enrich |> as.data.frame() |> rownames_to_column('tfs')
      
    if (enrichments_num == 1 ) {
      final_enrich <- matrix_enrich
    }
    else {
      final_enrich <- merge.data.frame(final_enrich , matrix_enrich , by = "tfs" , all = T)
    }
    }
    final_enrich_2 <- final_enrich |> column_to_rownames('tfs')


    saved_tfs <- c()

    for (column_name in colnames(final_enrich_2)) {
      col_df <- final_enrich_2[!is.na(final_enrich_2[[column_name]]), column_name, drop = FALSE]
      col_df <- col_df[order(-col_df[[column_name]]),, drop = FALSE]
      top_values <- head(col_df, top_n)
      saved_tfs <- c(saved_tfs, rownames(top_values))
    }

    saved_tfs <- unique(saved_tfs)

    # matrix for ordering
    cluster_map <- final_enrich_2[saved_tfs, ]
    cluster_map[is.na(cluster_map)] <- 0

    r.cluster <- fastcluster::hclust(dist(cluster_map, method = "euclidean"), method = "ward.D2")

    # -------------------------
    # 5. Long-format for plotting
    # -------------------------
    final_homer_scatterplot <- final_enrich_2[saved_tfs, ] |>
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
ggplot(final_homer_scatterplot, aes(x = Community, y = Factor)) +
      geom_point(aes(size = ComboScore),
                shape = 21, stroke = 0.6, color = "black") +
      geom_point(aes(size = ComboScore),
                shape = 21, fill = "black", stroke = 0.6, alpha = 0.4) +
      scale_x_discrete(drop=F)+
      theme_classic() +
      ylab("")+
  theme(text = element_text(size=10)) +
    theme(
      legend.position = "top",
      legend.text = element_text(size = 10)
    )
ggsave("/mnt/nas-safu02/sdigiove_workspace/check_th_TCHRONET/WGCNA/picture/monalisa_min100.pdf",device = cairo_pdf , 
  height = 18, 
  width = 10, 
  units = 'in',
  dpi = 300
)
