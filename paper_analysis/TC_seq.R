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

std_data <- tca@clusterRes@data

# 2. Define the clustering function for clusGap
# This function tells clusGap to use k-means
km_func <- function(x, k) {
  list(cluster = kmeans(x, centers = k, nstart = 25)$cluster)
}

# 3. Calculate the Gap Statistic
# K.max: maximum number of clusters to consider
# B: Number of Monte Carlo samples (use B=50 for a quick look, B=500 for publication)
gap_stat <- clusGap(std_data, 
                    FUNcluster = km_func, 
                    K.max = 8, 
                    B = 50, 
                    verbose = FALSE)



# 4. Visualize the results
fviz_gap_stat(gap_stat) +
  labs(title = "Gap Statistic for TCseq Clustering (k-means)")

tca <- timeclust(tca , algo = 'km' , k = 6 , standardize = TRUE)
timeclustplot(tca , col = 3)
save_files_community(tca , "/mnt/nas-safu02/sdigiove_workspace/check_th_TCHRONET/TCseq/km")

std_data <- tca@clusterRes@data
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
gap_stat_hc <- clusGap(std_data, 
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
