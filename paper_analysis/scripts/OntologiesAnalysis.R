suppressMessages (expr = { 
    library(rGREAT)
    library(tidyverse)
    library(clustree)
    library(GSEABase)
    library(clusterProfiler)
    library(enrichR)
    library(msigdbr)
})

# FUNCTIONS
interrogate_GREAT <- function(path_comm , comm_num , species , ontologies , th_qval = 0.05 , th_fc = 2 )  {
    gr = read_delim(path_comm,  col_names = c('chr','start','end') , delim = "\t" )#  

    # Create the GRanges
    gr = makeGRangesFromDataFrame(gr)

    # Define community name
    community_name = paste('community' , as.character(comm_num) , sep = "") 

    # Create job for great
    job = submitGreatJob(gr, version = "4.0.4" ,species = species , rule = "basalPlusExt" ,  adv_upstream = 5.0 , adv_downstream = 1.0 , request_interval = 12 ) #, adv_spa = 500

    # Obtain GREAT results
    tbl = getEnrichmentTables(job , ontology = ontologies )
    
    
    genes  = getRegionGeneAssociations(job)    
    target_genes = unname(unlist(as.data.frame(genes)$annotated_genes))

    GREAT_results = list()

    for (ontology in ontologies) {
        GREAT_results[[ontology]] = tbl[[ontology]]
    }
    

    columnNames <- paste(c("Binom_Adjp_BH" , "Binom_Fold_Enrichment") , community_name , sep="_")
    columnNames <- c('ID' , 'name' , columnNames)
    

    for (ontology_name in names(GREAT_results)) {
        GREAT_results[[ontology_name]]  <- GREAT_results[[ontology_name]] |>
        #dplyr::filter(Binom_Adjp_BH <= 0.01) |>
        dplyr::filter(Binom_Adjp_BH <= th_qval)|>
        dplyr::filter(Hyper_Adjp_BH <= th_qval) |>
        dplyr::filter(Binom_Fold_Enrichment >= th_fc) |>
        dplyr::select(ID , name , Binom_Adjp_BH , Binom_Fold_Enrichment)
        colnames( GREAT_results[[ontology_name]]) <- columnNames

    }
    return (list(GREAT_results , target_genes))
}

plotting_ontologies <- function(df , top_n = 10 , community_number)  {
  
  df_OnlyPval <- df %>% 
    dplyr::select( ID , name , paste('Binom_Adjp_BH_community' , community_number , sep = "")) %>% 
    gather(value = "Hyper_Adjp_BH" , key = "community" , -ID , -name)
  df_OnlyPval_TopValues <- df_OnlyPval %>% 
    group_by(community) %>% 
    arrange(Hyper_Adjp_BH) %>% 
    slice(1:top_n) 
  
  ggplot(df_OnlyPval_TopValues , aes(x = community , y = ID , color = -log10(Hyper_Adjp_BH))) +
    geom_point(size = 5) +
    theme_classic() +
    theme()
  
  clusteing_df = df_OnlyPval_TopValues %>% 
    tidyr::spread(key = community ,value = Hyper_Adjp_BH) %>% column_to_rownames("name") %>% dplyr::select(-ID)
  clusteing_df[is.na(clusteing_df)] = 0
  clusteing_df = clusteing_df[rowSums(clusteing_df != 0 ) >=1,]
  clusteing_df = clusteing_df 
  row.order <- stats::hclust(stats::dist(clusteing_df , method = 'euclidean'),method = "ward.D2")$order
  col.order <- stats::hclust(stats::dist(t(clusteing_df), method = 'euclidean'), method = "ward.D2")$order
  
  clusteing_df = clusteing_df[row.order, col.order]
  
  ggplot(data = df_OnlyPval_TopValues, aes(x = community, y = name)) + 
    geom_point(aes( color = -log10(Hyper_Adjp_BH)  ), size = 5  ) + #color =Percent  , , size = 5
    scale_color_gradient(low = "blue" , high = "red" , na.value ="white")+
    scale_radius() +#trans = "log2"
    theme_classic() +
    scale_y_discrete(limits=rownames(clusteing_df))+
    scale_x_discrete(limits =colnames(clusteing_df), guide = guide_axis(angle = 90))
}
GREAT_target_genes  <- function (file_list , ontologies_list) {
    final_target_gens = list()
    for (comm_num in  seq_along(list_files_communities)){
    file_path = list_files_communities[comm_num]
    # read the bed file
    
    GREAT_results = interrogate_GREAT(file_path ,  comm_num , "hg19" , ontologies)
    
    if (comm_num == 1 )  {
        Final_ontologies_dfs = GREAT_results[[1]]
        final_target_gens[[1]] = GREAT_results[[2]]
    }
    else {
        final_target_gens[[length(final_target_gens)+1]] = GREAT_results[[2]]
        for (ontology_name in names(GREAT_results[[1]])) {
        Final_ontologies_dfs[[ontology_name]]  <-  merge.data.frame(Final_ontologies_dfs[[ontology_name]]  , GREAT_results[[1]][[ontology_name]]  , by = c('ID' , 'name') , all = TRUE)
        }
    }
    }
    return( list(final_target_gens , Final_ontologies_dfs) )
}


# THP-1
ontologies <- c("GO Molecular Function","GO Biological Process","GO Cellular Component","Mouse Phenotype" ,"Mouse Phenotype Single KO","Human Phenotype")
list_files_communities = list.files(path = "/mnt/nas-safu/analysis/PhDsdigiove/method_coAcces/data/CellReport/communities_k27/res_1_5/beds/" , full.names = TRUE)

results  <- GREAT_target_genes(list_files_communities , ontologies)

msigdbr_collections(db_species = "Hs")
genesets <- msigdbr(species = "Homo sapiens" , db_species = 'HS', collection = "H") #, subcollection = "CP:REACTOME"
genesets_removed  <- genesets |> select(gs_name ,gene_symbol )
i = 1
j=0
for (x in reults[[1]]) {
    x_enrichr <- enricher( x , TERM2GENE = genesets_removed ) #   , universe = unname(unlist(reults[[1]]))
    x_df  <- x_enrichr@result |> filter(qvalue < 0.05) |> arrange ( -FoldEnrichment )

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

clusteing_df = final_df_selected_Columns |> select(Description , FoldEnrichment , community) |> spread(community , FoldEnrichment) |> column_to_rownames("Description")
clusteing_df[is.na(clusteing_df)] = 0
row.order <- stats::hclust(stats::dist(clusteing_df))$order
col.order <- stats::hclust(stats::dist(t(clusteing_df)))$order
clusteing_df[clusteing_df == 0] <-  NA
clusteing_df = clusteing_df[row.order, col.order]

final_df_selected_Columns  <- final_df_selected_Columns |> mutate(name = fct_relevel(Description, 
            rownames(clusteing_df)))
final_df_selected_Columns$community  <- factor(final_df_selected_Columns$community , levels = sort(final_df_selected_Columns |> pull(community) |> unique()))

ggplot(final_df_selected_Columns , aes(x = community , y = Description , color = qval_log10 , size = FoldEnrichment)) +
  geom_point() + #color =Percent  , , size = 5
  scale_color_gradient(low = "blue" , high = "red" , na.value ="white")+
  scale_radius() +#trans = "log2"
  theme_classic() +
  # scale_y_discrete(limits=rownames(clusteing_df)) +
  theme(legend.position="top" , legend.text = element_text(size=10)) 
  #scale_x_discrete(limits =colnames(clusteing_df), guide = guide_axis(angle = 90)) +
  # theme(axis.title.y=element_blank(), #axis.text.y = element_text(color = "grey20", size = 10, angle = 0, hjust = 1, vjust = 0, face = "plain") ,
  #       axis.text.y=element_blank(),
  #       legend.position="none")
#ggsave("/mnt/nas-safu/analysis/PhDsdigiove/method_coAcces/data/CellReport/pictures/HALLMARKS_withLabelswithLegends.png", width = 12, height = 5, units = "in", dpi = 300)


# LIVER DEVELOPMENT
ontologies <- c("GO Molecular Function","GO Biological Process","GO Cellular Component","Mouse Phenotype" ,"Mouse Phenotype Single KO","Human Phenotype")

list_files_communities = list.files(path = "/mnt/nas-safu/analysis/PhDsdigiove/method_coAcces/data/liver_devel_mouse_ENCODE/communities/res_0_8/beds/" , full.names = TRUE)
final_target_gens = list()

for (comm_num in  seq_along(list_files_communities)){
  file_path = list_files_communities[comm_num]
  print(file_path)
  # read the bed file
  
  GREAT_results = interrogate_GREAT(file_path ,  comm_num , "mm10" , ontologies)
  
  if (comm_num == 1 )  {
    Final_ontologies_dfs = GREAT_results[[1]]
    final_target_gens[[1]] = GREAT_results[[2]]
  }
  else {
    final_target_gens[[length(final_target_gens)+1]] = GREAT_results[[2]]
    for (ontology_name in names(GREAT_results[[1]])) {
      Final_ontologies_dfs[[ontology_name]]  <-  merge.data.frame(Final_ontologies_dfs[[ontology_name]]  , GREAT_results[[1]][[ontology_name]]  , by = c('ID' , 'name') , all = TRUE)
    }
  }
}

msigdbr_collections(db_species = "Mm")
genesets <- msigdbr(species = "Mus musculus" , db_species = 'MM', collection = "M8" ) # , subcollection = "GO:BP"
genesets_removed  <- genesets |> select(gs_name ,gene_symbol )

i = 1
for (x in final_target_gens) {
    x_enrichr <- enricher( x , TERM2GENE = genesets_removed , universe = unname(unlist(final_target_gens)) ) #  
    x_df  <- x_enrichr@result |> filter(qvalue < 0.05) |> arrange ( -FoldEnrichment )
    x_df  <- top_n(x_df,5,FoldEnrichment)
    x_df['community'] <- i

    if (i == 1) {
        final_df <- x_df
    } else {
        final_df <- rbind(final_df , x_df)
    }
    i = i + 1
}
clusteing_df = final_df_selected_Columns |> select(Description , FoldEnrichment , community) |> spread(community , FoldEnrichment) |> column_to_rownames("Description")

clusteing_df[is.na(clusteing_df)] = 0
row.order <- stats::hclust(stats::dist(clusteing_df))$order
clusteing_df[clusteing_df == 0] <-  NA
clusteing_df = clusteing_df[row.order, col.order]

final_df_selected_Columns  <- final_df_selected_Columns |> mutate(name = fct_relevel(Description, 
            rownames(clusteing_df)))

ggplot(final_df_selected_Columns , aes(x = community , y = Description , color = qval_log10 , size = FoldEnrichment)) +
  geom_point() + #color =Percent  , , size = 5
  scale_color_gradient(low = "blue" , high = "red" , na.value ="white")+
  scale_radius() +#trans = "log2"
  theme_classic() +
  scale_y_discrete(limits=rownames(clusteing_df)) +
  #theme(legend.position="top" , legend.text = element_text(size=10)) 
  scale_x_discrete(limits =colnames(clusteing_df), guide = guide_axis(angle = 90)) +
  theme(axis.title.y=element_blank(), #axis.text.y = element_text(color = "grey20", size = 10, angle = 0, hjust = 1, vjust = 0, face = "plain") ,
        axis.text.y=element_blank(),
        legend.position="none")
#ggsave("/mnt/nas-safu/analysis/PhDsdigiove/method_coAcces/data/CellReport/pictures/HALLMARKS_withLabelswithLegends.png", width = 12, height = 5, units = "in", dpi = 300)
# ggsave("/mnt/nas-safu/analysis/PhDsdigiove/method_coAcces/data/CellReport/pictures/HALLMARKS_NoLabelsNoLegends.png", width = 7, height = 5, units = "in", dpi = 300)

# B-ALL
ontologies <- c("GO Molecular Function","GO Biological Process","GO Cellular Component","Mouse Phenotype" ,"Mouse Phenotype Single KO","Human Phenotype")

list_files_communities = list.files(path = "/mnt/nas-safu/analysis/PhDsdigiove/method_coAcces/data/BALL/communities_broad/res_0_9/beds/" , full.names = TRUE)
final_target_gens = list()
for (comm_num in  seq_along(list_files_communities)){
  file_path = list_files_communities[comm_num]
  print(file_path)
  # read the bed file
  
  GREAT_results = interrogate_GREAT(file_path ,  comm_num , "hg19" , ontologies )
  
  if (comm_num == 1 )  {
    Final_ontologies_dfs = GREAT_results[[1]]
    final_target_gens[[1]] = GREAT_results[[2]]
  }
  else {
    final_target_gens[[length(final_target_gens)+1]] = GREAT_results[[2]]
    for (ontology_name in names(GREAT_results[[1]])) {
      Final_ontologies_dfs[[ontology_name]]  <-  merge.data.frame(Final_ontologies_dfs[[ontology_name]]  , GREAT_results[[1]][[ontology_name]]  , by = c('ID' , 'name') , all = TRUE)
    }
  }
}

msigdbr_collections(db_species = "HS")
genesets <- msigdbr(species = "Homo sapiens" , db_species = 'HS', collection = "C4" , subcollection = "3CA" ) #  , collection = "C2"   
genesets_removed  <- genesets |> select(gs_name ,gene_symbol )
i = 1
j=0
for (x in final_target_gens) {
    x_enrichr <- enricher( unique(x) , TERM2GENE = genesets_removed , universe = unique(unname(unlist(genesets_removed))) )  
    x_df  <- x_enrichr@result |> filter(qvalue < 0.1) #|> arrange ( qvalue )

    if(nrow(x_df) == 0) {    
        i = i + 1
        next
    }

    x_df  <- slice_min(x_df,FoldEnrichment, n = 5)
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
final_df_selected_Columns  <- final_df_selected_Columns |> mutate(name = fct_relevel(Description, 
            rownames(clusteing_df)))
final_df_selected_Columns$community  <- factor(final_df_selected_Columns$community , levels = sort(final_df_selected_Columns |> pull(community) |> unique()))

ggplot(final_df_selected_Columns , aes(x = community , y = Description , color = qval_log10 , size = FoldEnrichment)) +
  geom_point() + #color =Percent  , , size = 5
  scale_color_gradient(low = "blue" , high = "red" , na.value ="white")+
  scale_radius() +#trans = "log2"
  theme_classic() +
  # scale_y_discrete(limits=rownames(clusteing_df)) +
  #theme(legend.position="top" , legend.text = element_text(size=10)) 
  #scale_x_discrete(limits =colnames(clusteing_df), guide = guide_axis(angle = 90)) +
  theme(axis.title.y=element_blank(), #axis.text.y = element_text(color = "grey20", size = 10, angle = 0, hjust = 1, vjust = 0, face = "plain") ,
        axis.text.y=element_blank(),
        legend.position="none")
# ggsave("/mnt/nas-safu/analysis/PhDsdigiove/method_coAcces/data/BALL/communities_broad/res_0_9/msigdb_3CA.png", width = 12, height = 5, units = "in", dpi = 300)