
plotting_ontologies <- function( df , top_n = 10 )  {
  #### FUNCTION TO PLOT ONTOLOGIES FOR EACH COMMUNITY
  #### Parameters: 
  #### df : dataframe of Rgreat analysis
  #### top_n : specify how many top ontologies to be plotted for each community
  #### Return :
  #### ggplot dot plot of communities
  
  
  df_OnlyPval <- df %>% 
    gather(value = "Hyper_Adjp_BH" , key = "community" , -ID , -name)
  df_OnlyPval_TopValues <- df_OnlyPval %>% 
    group_by(community) %>% 
    arrange(Hyper_Adjp_BH) %>% 
    slice(1:top_n) 
  
  # ggplot(df_OnlyPval_TopValues , aes(x = community , y = ID , color = -log10(Hyper_Adjp_BH))) +
  #   geom_point(size = 5) +
  #   theme_classic() +
  #   theme()
  
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

interrogate_GREAT <- function(path_comm , comm_num , species , ontologies , th_qval = 0.01 , th_fc = 2 )  {
  #### FUNCTION TO INTERROGATE GREAT AND RETRIEVE THE ENRICHMENT FOR EACH ONTOLOGY
  #### Parameters: 
  #### path_comm : path of community's bed file
  #### comm_num : number of the community
  #### species : specify the species between 'hg' or 'mm'
  #### ontologies : vector with ontologies to be investigated through GREAT
  #### th_qval : qval threshold
  #### th_qval : fold change threshold
  #### Return :
  #### list of dataframes for each ontology
  
  
  gr = read_delim(path_comm,  col_names = c('chr','start','end') , delim = "\t" )#  
  
  # Create the GRanges
  gr = makeGRangesFromDataFrame(gr)
  
  # Define community name
  community_name = paste('community' , as.character(comm_num) , sep = "") 
  
  # Create job for great
  job = submitGreatJob(gr, version = "4.0.4" ,species = species , rule = "basalPlusExt" ,  adv_upstream = 2.0 , adv_downstream = 1.0 , request_interval = 12)
  
  # Obtain GREAT results
  tbl = getEnrichmentTables(job , ontology = ontologies )
  
  GREAT_results = list()
  
  for (ontology in ontologies) {
    GREAT_results[[ontology]] = tbl[[ontology]]
  }
  
  
  columnNames <- paste(c("Binom_Adjp_BH" , "Binom_Fold_Enrichment") , community_name , sep="_")
  columnNames <- c('ID' , 'name' , columnNames)
  
  
  for (ontology_name in names(GREAT_results)) {
    GREAT_results[[ontology_name]]  <- GREAT_results[[ontology_name]] |>
      dplyr::filter(Binom_Adjp_BH <= 0.01) |>
      dplyr::filter(Binom_Fold_Enrichment >= 2) |>
      dplyr::select(ID , name , Binom_Adjp_BH , Binom_Fold_Enrichment)
    colnames( GREAT_results[[ontology_name]]) <- columnNames
    
  }
  return (GREAT_results)
}

