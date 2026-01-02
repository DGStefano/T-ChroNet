#' Annotate regions according to a provided bed annotation file 
#'
#' @import dplyr
#' @import BSgenome
#' @import GenomicRanges
#' @import monaLisa
#' @import tidyr
#' 
#' @param object TChroNetNetwork object
#' @param negLog10Padj_th negative log 10 p-value threshold. defult = 4.0
#' @param log2enr_th log 2 enrichment threshold. defult = 1.0
#' @param top_n number of most enriched transcription factors for each community to be plotted
#' @return e ggplot2 object
#' @export
plot_enrichments_tfs <- function(object , negLog10Padj_th  = 4.0 , log2enr_th = 1.0 , top_n = 10) {
  if (!inherits(object, "TCrhoNetNetwork")) {
    stop("The 'object' must be a TCrhoNetNetwork instance.")
  }

  for (enrichments_num in seq_along(object@TFs)) {
    enrichments <- object@TFs[[enrichments_num]]
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
    p <- ggplot(final_homer_scatterplot, aes(x = Community, y = Factor)) +
      geom_point(aes(size = ComboScore),
                shape = 21, stroke = 0.6, color = "black") +
      geom_point(aes(size = ComboScore),
                shape = 21, fill = "black", stroke = 0.6, alpha = 0.4) +
      scale_x_discrete(drop=F)+
      theme_classic() +
      ylab("")
  
  return(p)

}