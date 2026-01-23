library(TCseq)
library(cluster)
library(factoextra)
library(dplyr)
library(tidyverse)
library(ggalluvial)
library(msigdbr)

get_great_genes <- function(cluster_df, coord_col, cluster_col, species = "mm10") {
  library(GenomicRanges)
  library(rGREAT)
  library(tidyverse)
  
  gene_list <- list()
  communities <- sort(unique(cluster_df[[cluster_col]]))
  
  for (comm in communities) {
    message(paste("Interrogating GREAT for cluster:", comm))
    
    # Format regions
    regions <- cluster_df[cluster_df[[cluster_col]] == comm, ] %>% 
      select(all_of(coord_col)) %>% 
      separate(!!sym(coord_col), c('chr','start','end'), sep ="-")
    
    gr <- makeGRangesFromDataFrame(regions)
    
    # Submit job
    job <- submitGreatJob(gr, version = "4.0.4", species = species, 
                          rule = "basalPlusExt", request_interval = 12)
    
    # Associate genes
    associations <- getRegionGeneAssociations(job)
    gene_list[[as.character(comm)]] <- unlist(unname(associations$annotated_genes))
  }
  return(gene_list)
}

run_pathway_enrichment <- function(target_genes_list, genesets_db, top_n_paths = 5) {
  library(clusterProfiler)
  library(DOSE)
  library(tidyverse)
  
  final_df <- data.frame()
  
  for (i in seq_along(target_genes_list)) {
    comm_name <- names(target_genes_list)[i]
    genes <- target_genes_list[[i]]
    
    if (length(genes) == 0) next
    
    # Run Enrichment
    x_enrichr <- enricher(genes, TERM2GENE = genesets_db)
    
    if (is.null(x_enrichr) || nrow(x_enrichr@result) == 0) next
    
    # Process results
    x_df <- x_enrichr@result %>% 
      filter(qvalue < 0.05) %>%
      mutate(
        FoldEnrichment = (parse_ratio(GeneRatio)) / (parse_ratio(BgRatio)),
        community = comm_name
      ) %>%
      arrange(-FoldEnrichment) %>%
      slice_max(FoldEnrichment, n = top_n_paths, with_ties = FALSE)
    
    final_df <- rbind(final_df, x_df)
  }
  return(final_df)
}

save_files_community <- function(tca, outdir ) {
  communities <- sort(unique(tca@clusterRes@cluster))

  for (community in communities) {
    peaks_list <- names(tca@clusterRes@cluster[tca@clusterRes@cluster == community])
    comm_df <- tca@genomicFeature |> filter(id %in% peaks_list) |> select(-id)
    outpath <- paste(outdir , "/community_" , as.character(community), ".bed" , sep = "")
    write_delim(comm_df , outpath , delim = "\t" , col_names = FALSE)
  }
}
species_genome <- "mm10"

dir <- ""
gf <- peakreference(dir = dir , patter = "bed")

bamfiles <- data.frame(
    list(
        sampleid = c('t11_rep1','t11_rep2','t12_rep1','t12_rep2','t13_rep1','t13_rep2','t14_rep1','t14_rep2','t15_rep1','t15_rep2','t16_rep1','t16_rep2'),
        timepoint = c('t11','t11','t12','t12','t13','t13','t14','t14','t15','t15','t16','t16'),
        group = c("1","1","2","2","3","3","4","4","5","5","6","6"),
        bamfile = c("liver_mouse_t11_rep1.bam", "liver_mouse_t11_rep2.bam","liver_mouse_t12_rep1.bam", "liver_mouse_t12_rep2.bam","liver_mouse_t13_rep1.bam", "liver_mouse_t13_rep2.bam","liver_mouse_t14_rep1.bam", "liver_mouse_t14_rep2.bam","liver_mouse_t15_rep1.bam", "liver_mouse_t15_rep2.bam","liver_mouse_t16_rep1.bam", "liver_mouse_t16_rep2.bam")
    )
)

tca <- TCA(design = bamfiles , genomicFeature = gf)

dir.BAM <- ""
tca <- countReads(tca , dir = dir.BAM)

tca <- DBanalysis(tca)

tca <- timecourseTable(tca)

t <- tcTable(tca)


# Find best cluster through silhouette score
res_k <- fviz_nbclust(t, kmeans, method = "silhouette", k.max = 5)

# Extract the k with the highest average silhouette width
best_k <- as.numeric(res_k$data$clusters[which.max(res_k$data$y)])

# Final clustering with the selected best_k
tca <- timeclust(tca, algo = 'km', k = best_k, standardize = TRUE)
timeclustplot(tca , col = 3)
save_files_community(tca , "~/T-ChroNet/paper_analysis/data/banchmark/tcseq_km")


t <- tcTable(tca)
# std_data should be your standardized matrix (Z-scores)
# method = "silhouette"
# FUNcluster = hcut (this performs hierarchical clustering)
res_k <- fviz_nbclust(t, 
                      FUNcluster = hcut, 
                      method = "silhouette", 
                      hc_method = "ward.D2", # You can choose ward.D2, complete, etc.
                      k.max = 10)

# Retrieve the optimal K
best_k <- as.numeric(res_k$data$clusters[which.max(res_k$data$y)])

# Final clustering with the selected best_k
tca <- timeclust(tca, algo = 'hc', k = best_k, standardize = TRUE)
timeclustplot(tca , col = 1)
save_files_community(tca , "~/T-ChroNet/paper_analysis/data/banchmark/tcseq_hc")

#### COMPARISON HC
files_lists <- list.files("~/T-ChroNet/paper_analysis/data/banchmark/tcseq_hc" , full.names = T)
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


files_lists <- list.files("~/T-ChroNet/paper_analysis/data/LiverDevelopment/communities/" , full.names = T)
i = 0
for ( x in  seq_along(files_lists)) {
  df <- read_delim(files_lists[[x]] , delim = "\t" , col_names = F) |> unite('peaks' , c(X1,X2,X3) , sep ="-")
  df['cluster_TCNR'] <- x
  if ( i == 0 ){
    final_tcnr_clst <- df
    i = 1 
  }
  else {
    final_tcnr_clst <- rbind(final_tcnr_clst , df)
  }
}

#### COMPARISON KM
files_lists <- list.files("~/T-ChroNet/paper_analysis/data/banchmark/tcseq_km/" , full.names = T)
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

final_comparison <- merge.data.frame(final_tcseq_clst_hc , final_tcnr_clst , by = 'peaks')
final_comparison <- merge.data.frame(final_comparison , final_tcseq_clst_km , by = 'peaks')

final_comparison[is.na(final_comparison)] <- 100
final_comparison$cluster_Tcseq_hc <- factor(final_comparison$cluster_Tcseq_hc , levels = sort(unique(final_comparison$cluster_Tcseq_hc)))
final_comparison$cluster_TCNR <- factor(final_comparison$cluster_TCNR , levels = sort(unique(final_comparison$cluster_TCNR)))
final_comparison$cluster_Tcseq_km <- factor(final_comparison$cluster_Tcseq_km , levels = sort(unique(final_comparison$cluster_Tcseq_km)))


pdf("~/T-ChroNet/paper_analysis/data/banchmark/TC_seq_hc_our_km_LIVER.pdf"  ,
  height = 5, 
  width = 9,
  family = "ArialMT",
  useDingbats = FALSE)
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
dev.off()


##### ENRICHMENTS
genes_km  <- get_great_genes(final_tcseq_clst_km, "peaks", "cluster_Tcseq_km", species_genome)
genes_hc  <- get_great_genes(final_tcseq_clst_hc, "peaks", "cluster_Tcseq_hc", species_genome)
genes_tcn <- get_great_genes(final_tcnr_clst,     "peaks", "cluster_TCNR",     species_genome)

# Apply enrichment function
genesets <- msigdbr(species = "Mus musculus" , db_species = 'MM', collection = "M8" )
genesets_removed  <- genesets |> dplyr::select(gs_name ,gene_symbol )


final_df_km  <- run_pathway_enrichment(genes_km,  genesets_removed)
final_df_hc  <- run_pathway_enrichment(genes_hc,  genesets_removed)
final_df_tcn <- run_pathway_enrichment(genes_tcn, genesets_removed)

#######
# 1. Add a source column to each dataframe and combine them
final_df_km$source <- "K-means"
final_df_hc$source <- "Hierarchical"
final_df_tcn$source <- "TCN"

combined_df <- rbind(final_df_km, final_df_hc, final_df_tcn)

# 2. Select columns and calculate -log10 q-value
final_df_selected_Columns <- combined_df %>% 
  select(source, community, Description, FoldEnrichment, qvalue) %>%
  mutate(qval_log10 = -log10(qvalue))

# 3. Clean Description strings (removing prefixes and underscores)
final_df_selected_Columns$Description <- final_df_selected_Columns$Description %>%
  str_remove("DESCARTES_") %>%
  str_remove("TABULA_MURIS") %>%
  str_replace_all("_", " ")

# 4. Clustering/Ordering Logic (Global ordering for the Y-axis)
# This ensures that similar terms across all datasets are grouped together
clustering_df <- final_df_selected_Columns %>% 
  select(Description, FoldEnrichment, community, source) %>%
  # Create a unique ID for community within each source to avoid column collisions
  mutate(unique_comm = paste0(source, "_", community)) %>%
  select(-community, -source) %>%
  pivot_wider(names_from = unique_comm, values_from = FoldEnrichment, values_fill = 0) %>%
  column_to_rownames("Description")

row_order <- stats::hclust(stats::dist(clustering_df))$order
ordered_terms <- rownames(clustering_df)[row_order]

# Apply the order and factor levels
final_df_selected_Columns <- final_df_selected_Columns %>%
  mutate(
    Description = factor(Description, levels = ordered_terms),
    source = factor(source, levels = c("K-means", "Hierarchical", "TCN")), # Order of facets
    community = factor(community)
  )

pdf("~/T-ChroNet/paper_analysis/data/banchmark/CellAnnotations_alldf_TCseq.png"  ,
  height = 5, 
  width = 9,
  family = "ArialMT",
  useDingbats = FALSE)
# 5. Plotting with Facets
ggplot(final_df_selected_Columns, aes(x = community, y = Description, color = qval_log10, size = FoldEnrichment)) +
  geom_point() +
  scale_color_gradient(low = "blue", high = "red", na.value = "white") +
  scale_radius() +
  theme_bw() + 
  theme(
    strip.background = element_blank(),
    legend.position = "top",
    legend.text = element_text(size = 10),
    text = element_text(size = 10),
    panel.spacing = unit(1, "lines"), # This controls the "space" between datasets
    axis.text.x = element_text(angle = 0)
  ) +
  # Use facet_grid to separate datasets. 
  # scales = "free_x" and space = "free_x" are CRUCIAL to keep community widths equal
  facet_grid(. ~ source, scales = "free_x", space = "free_x") +
  xlab("Community") +
  ylab("") +
  labs(color = "-log10(q-value)", size = "Fold Enrichment")
dev.off()
#######

library(WGCNA)
library(tidyverse)
library(flashClust)

allowWGCNAThreads(nThreads = 40)
powers = c(c(1:10), seq(from = 12, to = 40, by = 2))

input_mat <- read_tsv("~/T-ChroNet/paper_analysis/data/LiverDevelopment/multicov_all_sites_all_timepoint.tsv") |> unite('peaks' , c(chromosome, start,end) , sep ="-") |> tibble::column_to_rownames('peaks')
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
dynamicMods_1000 = cutreeDynamic(dendro = geneTree,  method="tree", minClusterSize = 1000);
peak_clusters_1000 <- data.frame(list(peaks = colnames(input_mat) , clusuters_wgcna = dynamicMods_1000))
dynamicMods_100 = cutreeDynamic(dendro = geneTree,  method="tree", minClusterSize = 100);
peak_clusters_100 <- data.frame(list(peaks = colnames(input_mat) , clusuters_wgcna_100 = dynamicMods_100))
peak_clusters_100$peaks <- gsub("_", "-", peak_clusters_100$peaks)
peak_clusters_1000$peaks <- gsub("_", "-", peak_clusters_1000$peaks)

write_delim(peak_clusters_100 , "~/T-ChroNet/paper_analysis/data/banchmark/peaks_clust_100.tsv" , delim = "\t")
write_delim(peak_clusters_1000 , "~/T-ChroNet/paper_analysis/data/banchmark/peaks_clust_1000.tsv" , delim = "\t")

peak_clusters_100 <- read_delim("~/T-ChroNet/paper_analysis/data/banchmark/peaks_clust_100.tsv" , delim = "\t")
peak_clusters_1000 <- read_delim("~/T-ChroNet/paper_analysis/data/banchmark/peaks_clust_1000.tsv" , delim = "\t")


final_comparison <- merge.data.frame(peak_clusters_100 , final_tcnr_clst , by = 'peaks')
final_comparison <- merge.data.frame(final_comparison , peak_clusters_1000 , by = 'peaks')

final_comparison$clusuters_wgcna <- factor(final_comparison$clusuters_wgcna , levels = sort(unique(final_comparison$clusuters_wgcna)))
final_comparison$cluster_TCNR <- factor(final_comparison$cluster_TCNR , levels = sort(unique(final_comparison$cluster_TCNR)))
final_comparison$clusuters_wgcna_100 <- factor(final_comparison$clusuters_wgcna_100 , levels = sort(unique(final_comparison$clusuters_wgcna_100)))

pdf("~/T-ChroNet/paper_analysis/data/banchmark/wgcna1000_hc_our_wgcna100_LIVER.pdf"  ,
  height = 5, 
  width = 9,
  family = "ArialMT",
  useDingbats = FALSE)
# Plotting the transition between 'clusters' and 'clusters_tcnr'
ggplot(final_comparison, aes(axis1 = clusuters_wgcna, axis2 = cluster_TCNR , axis3 = clusuters_wgcna_100)) +
  geom_alluvium(aes(fill = cluster_TCNR), width = 1/12) +
  geom_stratum(width = 1/12, fill = "grey80", color = "black") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("TCseq Clusters HC", "TCNR Clusters" ,"TCseq Clusters KM"), expand = c(.05, .05 , .05)) +
  labs(title = "Peaks Cluster Membership Transition",
       y = "Number of Peaks") +
  theme_minimal() +
  theme(text = element_text(size = 15))
dev.off()



genes_wgcna <- get_great_genes(peak_clusters_1000, "peaks", "clusuters_wgcna", species_genome)
final_df_wgcna <- run_pathway_enrichment(genes_wgcna, genesets_removed)
final_df_wgcna$source <- "WGCNA"


combined_df <- rbind(final_df_wgcna, final_df_tcn)

# 2. Select columns and calculate -log10 q-value
final_df_selected_Columns <- combined_df %>% 
  select(source, community, Description, FoldEnrichment, qvalue) %>%
  mutate(qval_log10 = -log10(qvalue))

# 3. Clean Description strings (removing prefixes and underscores)
final_df_selected_Columns$Description <- final_df_selected_Columns$Description %>%
  str_remove("DESCARTES_") %>%
  str_remove("TABULA_MURIS") %>%
  str_replace_all("_", " ")

# 4. Clustering/Ordering Logic (Global ordering for the Y-axis)
# This ensures that similar terms across all datasets are grouped together
clustering_df <- final_df_selected_Columns %>% 
  select(Description, FoldEnrichment, community, source) %>%
  # Create a unique ID for community within each source to avoid column collisions
  mutate(unique_comm = paste0(source, "_", community)) %>%
  select(-community, -source) %>%
  pivot_wider(names_from = unique_comm, values_from = FoldEnrichment, values_fill = 0) %>%
  column_to_rownames("Description")

row_order <- stats::hclust(stats::dist(clustering_df))$order
ordered_terms <- rownames(clustering_df)[row_order]

# Apply the order and factor levels
final_df_selected_Columns <- final_df_selected_Columns %>%
  mutate(
    Description = factor(Description, levels = ordered_terms),
    source = factor(source, levels = c("K-means", "Hierarchical", "TCN")), # Order of facets
    community = factor(community)
  )

pdf("~/T-ChroNet/paper_analysis/data/banchmark/CellAnnotations_alldf_WGCNA.png"  ,
  height = 5, 
  width = 9,
  family = "ArialMT",
  useDingbats = FALSE)
# 5. Plotting with Facets
ggplot(final_df_selected_Columns, aes(x = community, y = Description, color = qval_log10, size = FoldEnrichment)) +
  geom_point() +
  scale_color_gradient(low = "blue", high = "red", na.value = "white") +
  scale_radius() +
  theme_bw() + 
  theme(
    strip.background = element_blank(),
    legend.position = "top",
    legend.text = element_text(size = 10),
    text = element_text(size = 10),
    panel.spacing = unit(1, "lines"), # This controls the "space" between datasets
    axis.text.x = element_text(angle = 0)
  ) +
  # Use facet_grid to separate datasets. 
  # scales = "free_x" and space = "free_x" are CRUCIAL to keep community widths equal
  facet_grid(. ~ source, scales = "free_x", space = "free_x") +
  xlab("Community") +
  ylab("") +
  labs(color = "-log10(q-value)", size = "Fold Enrichment")
dev.off()

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
                      matrixtype = "PWM",
                      all_versions = FALSE))

# disconnect Db
RSQLite::dbDisconnect(JASPARConnect)

genome_obj <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10

process_monalisa_enrichment <- function(cluster_df, 
                                        coord_col, 
                                        cluster_col, 
                                        pwms, 
                                        genome_obj, 
                                        method_label = "Method",
                                        negLog10Padj_th = 4.0, 
                                        log2enr_th = 1.0, 
                                        top_n = 20) {
  
  library(tidyverse)
  library(monaLisa)
  library(GenomicRanges)
  library(SummarizedExperiment)
  
  enrich_list <- list()
  clusters <- sort(unique(cluster_df[[cluster_col]]))
  
  for (i in clusters) {
    message(paste("Processing", method_label, "cluster:", i))
    
    nodes <- cluster_df[cluster_df[[cluster_col]] == i, ] %>%
      dplyr::select(all_of(coord_col)) %>%
      tidyr::separate(!!sym(coord_col), c('chr', 'start', 'end'), sep = "-")
    
    gr <- GenomicRanges::makeGRangesFromDataFrame(nodes)
    seqs <- BSgenome::getSeq(genome_obj, gr)
    
    se_genome <- monaLisa::calcBinnedMotifEnrR(
      seqs = seqs, pwmL = pwms, background = "genome",            
      genome = genome_obj, genome.oversample = 2,
      BPPARAM = BiocParallel::MulticoreParam(20)
    )
    enrich_list[[as.character(i)]] <- se_genome
  }
  
  final_enrich <- NULL
  
  for (idx in seq_along(enrich_list)) {
    enr <- enrich_list[[idx]]
    comm_name <- names(enrich_list)[idx]
    
    sel_padj <- apply(assay(enr, "negLog10Padj"), 1, function(x) max(abs(x), 0, na.rm = TRUE)) > negLog10Padj_th
    sel_enr  <- apply(assay(enr, "log2enr"), 1, function(x) max(x, na.rm = TRUE)) > log2enr_th
    
    seSel <- enr[sel_padj & sel_enr, ]
    
    mat <- as.data.frame(assay(seSel, "log2enr"))
    motif_meta <- as.data.frame(mcols(seSel)[c("motif.id", "motif.name")])
    
    # FIX: Create a unique ID to avoid the 'duplicate row.names' error
    # Format: "Arnt (MA0004.1)"
    rownames(mat) <- paste0(motif_meta$motif.name, " (", motif_meta$motif.id, ")")
    
    colnames(mat) <- paste0("community_", comm_name)
    mat <- mat %>% rownames_to_column("motif_unique")
    
    if (is.null(final_enrich)) {
      final_enrich <- mat
    } else {
      final_enrich <- full_join(final_enrich, mat, by = "motif_unique")
    }
  }
  
  final_matrix <- final_enrich %>% column_to_rownames("motif_unique")
  
  saved_motifs <- map(colnames(final_matrix), function(col) {
    final_matrix %>%
      filter(!is.na(.data[[col]])) %>%
      arrange(desc(.data[[col]])) %>%
      head(top_n) %>%
      rownames()
  }) %>% unlist() %>% unique()
  
  cluster_map <- final_matrix[saved_motifs, ]
  cluster_map[is.na(cluster_map)] <- 0
  r_cluster <- fastcluster::hclust(dist(cluster_map, method = "euclidean"), method = "ward.D2")
  
  long_df <- final_matrix[saved_motifs, ] %>%
    rownames_to_column("motif_unique") %>%
    pivot_longer(cols = -motif_unique, names_to = "Community", values_to = "ComboScore") %>%
    mutate(
      Community = as.numeric(gsub("\\D", "", Community)),
      Factor_Full = factor(motif_unique, levels = r_cluster$labels[r_cluster$order]),
      # Create a clean label for the plot (remove the ID)
      Factor = gsub(" \\(MA.*\\)", "", Factor_Full),
      Method = method_label
    )
  
  return(long_df)
}

# 1. Process all three
res_hc  <- process_monalisa_enrichment(final_tcseq_clst_hc, "peaks", "cluster_Tcseq_hc", pwms, genome_obj, "Hierarchical")
res_km  <- process_monalisa_enrichment(final_tcseq_clst_km, "peaks", "cluster_Tcseq_km", pwms, genome_obj, "K-means")
res_tcn <- process_monalisa_enrichment(final_tcnr_clst, "peaks", "cluster_TCNR", pwms, genome_obj, "TCN")

# 2. Combine
all_results <- bind_rows(res_hc, res_km, res_tcn)

# 3. Plot together with spacing
pdf("~/T-ChroNet/paper_analysis/data/banchmark/TFs_alldf_TCseq.pdf"  ,
  height = 8, 
  width = 7,
  family = "ArialMT",
  useDingbats = FALSE) 
ggplot(all_results, aes(x = factor(Community), y = Factor_Full)) +
  geom_point(aes(size = ComboScore), shape = 21, stroke = 0.6, color = "black") +
  geom_point(aes(size = ComboScore), shape = 21, fill = "black", stroke = 0.6, alpha = 0.4) +
  facet_grid(. ~ Method, scales = "free_x", space = "free_x") +
  scale_y_discrete(labels = function(x) gsub(" \\(MA.*\\)", "", x)) + # Cleans labels on Y axis
  theme_bw() +
  theme(
    panel.spacing = unit(2, "lines"), 
    axis.text.y = element_text(size = 7)
  ) +
  labs(x = "Community ID", y = "Transcription Factor", size = "Log2 Enr")
dev.off()
res_wgcna <- process_monalisa_enrichment(peak_clusters_1000, "peaks", "clusuters_wgcna", pwms, genome_obj, "WGCNA")

all_results <- bind_rows(res_wgcna, res_tcn)

# 3. Plot together with spacing
pdf("~/T-ChroNet/paper_analysis/data/banchmark/TFs_alldf_WGCNA.pdf"  ,
  height = 8, 
  width = 7,
  family = "ArialMT",
  useDingbats = FALSE)
ggplot(all_results, aes(x = factor(Community), y = Factor_Full)) +
  geom_point(aes(size = ComboScore), shape = 21, stroke = 0.6, color = "black") +
  geom_point(aes(size = ComboScore), shape = 21, fill = "black", stroke = 0.6, alpha = 0.4) +
  facet_grid(. ~ Method, scales = "free_x", space = "free_x") +
  scale_y_discrete(labels = function(x) gsub(" \\(MA.*\\)", "", x)) + # Cleans labels on Y axis
  theme_bw() +
  theme(
    panel.spacing = unit(2, "lines"), 
    axis.text.y = element_text(size = 7)
  ) +
  labs(x = "Community ID", y = "Transcription Factor", size = "Log2 Enr")
dev.off()
