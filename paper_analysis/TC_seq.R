library(TCseq)
library(cluster)
library(factoextra)
library(dplyr)
library(tidyverse)

save_files_community <- function(tca, outdir ) {
  communities <- sort(unique(tca@clusterRes@cluster))

  for (community in communities) {
    peaks_list <- names(tca@clusterRes@cluster[tca@clusterRes@cluster == community])
    comm_df <- tca@genomicFeature |> filter(id %in% peaks_list) |> select(-id)
    outpath <- paste(outdir , "/community_" , as.character(community), ".bed" , sep = "")
    write_delim(comm_df , outpath , delim = "\t" , col_names = FALSE)
  }
}


dir <- "/mnt/nas-safu01/analysis/PhDsdigiove/method_coAcces/data/liver_devel_mouse_ENCODE/bed/"
gf <- peakreference(dir = dir , patter = "bed")

bamfiles <- data.frame(
    list(
        sampleid = c('t11_rep1','t11_rep2','t12_rep1','t12_rep2','t13_rep1','t13_rep2','t14_rep1','t14_rep2','t15_rep1','t15_rep2','t16_rep1','t16_rep2'),
        timepoint = c('t11','t11','t12','t12','t13','t13','t14','t14','t15','t15','t16','t16'),
        group = c(1,1,2,2,3,3,4,4,5,5,6,6),
        bamfile = c("liver_mouse_t11_rep1.bam", "liver_mouse_t11_rep2.bam","liver_mouse_t12_rep1.bam", "liver_mouse_t12_rep2.bam","liver_mouse_t13_rep1.bam", "liver_mouse_t13_rep2.bam","liver_mouse_t14_rep1.bam", "liver_mouse_t14_rep2.bam","liver_mouse_t15_rep1.bam", "liver_mouse_t15_rep2.bam","liver_mouse_t16_rep1.bam", "liver_mouse_t16_rep2.bam")
    )
)

tca <- TCA(design = bamfiles , genomicFeature = gf)

dir.BAM <- "/mnt/nas-safu01/analysis/PhDsdigiove/method_coAcces/data/liver_devel_mouse_ENCODE/bam/"
tca <- countReads(tca , dir = dir.BAM)

tca <- DBanalysis(tca , filter.type = 'raw' , filter.value = 5  )

tca <- timecourseTable(tca , value = "expression" , control.group = "1" , norm.method = 'cpm' , filter = TRUE)

t <- tcTable(tca)
# 2. Define the clustering function for clusGap
# This function tells clusGap to use k-means
km_func <- function(x, k) {
  list(cluster = kmeans(x, centers = k, nstart = 1)$cluster)
}

# 3. Calculate the Gap Statistic
# K.max: maximum number of clusters to consider
# B: Number of Monte Carlo samples (use B=50 for a quick look, B=500 for publication)
gap_stat <- clusGap(t, 
                    FUNcluster = km_func, 
                    K.max = 8, 
                    B = 50, 
                    verbose = FALSE)

# method = "silhouette"
res_k <- fviz_nbclust(t, kmeans, method = "silhouette", k.max = 10)

# Extract the k with the highest average silhouette width
best_k <- as.numeric(res_k$data$clusters[which.max(res_k$data$y)])

# 4. Visualize the results
fviz_gap_stat(gap_stat) +
  labs(title = "Gap Statistic for TCseq Clustering (k-means)")

tca <- timeclust(tca , algo = 'km' , k = 6 , standardize = TRUE)
timeclustplot(tca , col = 3)
save_files_community(tca , "/mnt/nas-safu02/sdigiove_workspace/check_th_TCHRONET/TCseq/km")


t <- tcTable(tca)
# 1. Define the Hierarchical Clustering wrapper
hc_func <- function(x, k) {
  # Perform hierarchical clustering
  d <- dist(x)
  hc <- hclust(d, method = "complete") # Default TCseq linkage
  # Return cluster assignments
  list(cluster = cutree(hc, k = k))
}

# 2. Calculate the Gap Statistic
# B=50 is standard for testing; increase to 500 for final results
gap_stat_hc <- clusGap(t, 
                       FUNcluster = hc_func, 
                       K.max = 10, 
                       B = 50, 
                       verbose = FALSE)

# 3. Visualize
fviz_gap_stat(gap_stat_hc) +
  labs(title = "Gap Statistic: TCseq Hierarchical Clustering")

tca <- timeclust(tca , algo = 'hc' , k = 3 , standardize = TRUE)
timeclustplot(tca , col = 1)
save_files_community(tca , "/mnt/nas-safu02/sdigiove_workspace/check_th_TCHRONET/TCseq/hc")

#### COMPARISON HC
files_lists <- list.files("/mnt/nas-safu02/sdigiove_workspace/check_th_TCHRONET/TCseq/hc/" , full.names = T)
i = 0
for ( x in  seq_along(files_lists)) {
  df <- read_delim(files_lists[[x]] , delim = "\t" , col_names = F)
  df <- df |> unite('peaks' , c('X1','X2','X3') , sep ="-")
  df['cluster_Tcseq_hc'] <- x
  if ( i == 0 ){
    final_tcseq_clst_hc <- df
    i = 1 
  }
  else {
    final_tcseq_clst_hc <- rbind(final_tcseq_clst_hc , df)
  }
}


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

final_comparison <- merge.data.frame(final_tcseq_clst , final_tcnr_clst , by.x = 'peaks' , by.y = 'X1')

# Plotting the transition between 'clusters' and 'clusters_tcnr'
ggplot(final_comparison, aes(axis1 = cluster_Tcseq, axis2 = cluster_TCNR)) +
  geom_alluvium(aes(fill = cluster_Tcseq), width = 1/12) +
  geom_stratum(width = 1/12, fill = "grey80", color = "black") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("TCseq Clusters", "TCNR Clusters"), expand = c(.05, .05)) +
  labs(title = "Gene Cluster Membership Transition",
       y = "Number of Genes/Peaks") +
  theme_minimal()

#### COMPARISON KM
files_lists <- list.files("/mnt/nas-safu02/sdigiove_workspace/check_th_TCHRONET/TCseq/km/" , full.names = T)
i = 0
for ( x in  seq_along(files_lists)) {
  df <- read_delim(files_lists[[x]] , delim = "\t" , col_names = F)
  df <- df |> unite('peaks' , c('X1','X2','X3') , sep ="-")
  df['cluster_Tcseq_km'] <- x
  if ( i == 0 ){
    final_tcseq_clst_km <- df
    i = 1 
  }
  else {
    final_tcseq_clst_km <- rbind(final_tcseq_clst_km , df)
  }
}


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

final_comparison <- merge.data.frame(final_tcseq_clst_hc , final_tcnr_clst , by.x = 'peaks' , by.y = 'X1')
final_comparison <- merge.data.frame(final_comparison , final_tcseq_clst_km , by = 'peaks')

final_comparison$cluster_Tcseq_hc <- factor(final_comparison$cluster_Tcseq_hc , levels = sort(unique(final_comparison$cluster_Tcseq_hc)))
final_comparison$cluster_TCNR <- factor(final_comparison$cluster_TCNR , levels = sort(unique(final_comparison$cluster_TCNR)))
final_comparison$cluster_Tcseq_km <- factor(final_comparison$cluster_Tcseq_km , levels = sort(unique(final_comparison$cluster_Tcseq_km)))

# Plotting the transition between 'clusters' and 'clusters_tcnr'
ggplot(final_comparison, aes(axis1 = cluster_Tcseq_hc, axis2 = cluster_TCNR , axis3 = cluster_Tcseq_km)) +
  geom_alluvium(aes(fill = cluster_TCNR), width = 1/12) +
  geom_stratum(width = 1/12, fill = "grey80", color = "black") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("TCseq Clusters HC", "TCNR Clusters" ,"TCseq Clusters KM"), expand = c(.05, .05 , .05)) +
  labs(title = "Peaks Cluster Membership Transition",
       y = "Number of Peaks") +
  theme_minimal() +
  theme(text = element_text(size = 15))
ggsave("/mnt/nas-safu02/sdigiove_workspace/check_th_TCHRONET/TCseq/pictures/TC_seq_hc_our_km_LIVER.pdf" , height = 5 , width = 9 , units = 'in' , dpi = 300)


library(rGREAT)
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
for (comm in sort(unique(final_tcseq_clst_km$cluster_Tcseq_km))) {
  
  regions <- final_tcseq_clst_km[final_tcseq_clst_km$cluster_Tcseq_km == comm, ] |> select(peaks) |> separate(peaks,c('chr','start','end') , sep ="-")
    
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
    
  
  final_target_genes[[comm]] <- unlist(unname(GREAT_results$annotated_genes))
}

final_target_genes <- list()  
for (comm in sort(unique(final_tcseq_clst_hc$cluster_Tcseq_hc))) {
  
  regions <- final_tcseq_clst_hc[final_tcseq_clst_hc$cluster_Tcseq_hc == comm, ] |> select(peaks) |> separate(peaks,c('chr','start','end') , sep ="-")
    
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
    
  
  final_target_genes[[comm]] <- unlist(unname(GREAT_results$annotated_genes))
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
ggsave("/mnt/nas-safu02/sdigiove_workspace/check_th_TCHRONET/TCseq/pictures/Cell_annotations_TCseq_hc.pdf",device = cairo_pdf , 
  height = 5, 
  width = 9, 
  units = 'in',
  dpi = 300
)

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
  for (i in sort(unique(final_tcseq_clst_hc$cluster_Tcseq_hc))) {
    community_nodes <- final_tcseq_clst_hc[final_tcseq_clst_hc$cluster_Tcseq_hc == i, ] |> dplyr::select(peaks) |> tidyr::separate(peaks,c('chr','start','end') , sep ="-")

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

    r.cluster <- hclust(dist(cluster_map, method = "euclidean"), method = "ward.D2")

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
ggsave("/mnt/nas-safu02/sdigiove_workspace/check_th_TCHRONET/TCseq/pictures/monalisa_HC.pdf",device = cairo_pdf , 
  height = 9, 
  width = 5, 
  units = 'in',
  dpi = 300
)


monaLisa_enrichments <- list()
  for (i in sort(unique(final_tcseq_clst_km$cluster_Tcseq_km))) {
    community_nodes <- final_tcseq_clst_km[final_tcseq_clst_km$cluster_Tcseq_km == i, ] |> dplyr::select(peaks) |> tidyr::separate(peaks,c('chr','start','end') , sep ="-")

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

    r.cluster <- hclust(dist(cluster_map, method = "euclidean"), method = "ward.D2")

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
ggsave("/mnt/nas-safu02/sdigiove_workspace/check_th_TCHRONET/TCseq/pictures/monalisa_KM.pdf",device = cairo_pdf , 
  height = 9, 
  width = 5, 
  units = 'in',
  dpi = 300
)



monaLisa_enrichments <- list()
  for (i in sort(unique(final_tcnr_clst$cluster_TCNR))) {
    community_nodes <- final_tcnr_clst[final_tcnr_clst$cluster_TCNR == i, ] |> dplyr::select(X1) |> tidyr::separate(X1,c('chr','start','end') , sep ="-")

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
ggsave("/mnt/nas-safu02/sdigiove_workspace/check_th_TCHRONET/TCseq/pictures/monalisa_Original.pdf",device = cairo_pdf , 
  height = 9, 
  width = 5, 
  units = 'in',
  dpi = 300
)
