library(GenomicRanges)

library(BSgenome.Mmusculus.UCSC.mm10)
library(phastCons100way.UCSC.hg38)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(dplyr)
library(purrr)
library(tidyverse)
library(rstatix)

profile_community_biology_mm10 <- function(df_total) {
  
  # 1. Prepare Annotations
  genome <- BSgenome.Mmusculus.UCSC.mm10
  phast <- getGScores("phastCons60way.UCSC.mm10")
  txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
  all_tss <- promo_ranges <- promoters(txdb, upstream=0, downstream=1)
  returning_df <- data.frame()
  for ( community_name in unique(df_total$community)){
    df <- df_total |> filter(community == community_name)
    # Convert df to GRanges
    gr <- makeGRangesFromDataFrame(df, keep.extra.columns = TRUE , na.rm = TRUE)
    
    # --- Feature: Region Length ---
    df$region_length <- width(gr)
    
    # --- Feature: GC Content ---
    seqs <- getSeq(genome, gr)
    df$gc_content <- as.numeric(letterFrequency(seqs, letters = "GC", as.prob = TRUE))
    
    # --- Feature: TSS Distance ---
    # Find distance to the nearest Transcription Start Site
    dist_to_tss <- distanceToNearest(gr, all_tss)
    df$tss_distance <- mcols(dist_to_tss)$distance
    
    df$score_phylo <- score(phast, gr)
    # Keep only the new features and the community ID
    df$community_id <- community_name
    returning_df <- rbind(returning_df , df)
  }
  return(returning_df %>% select( chromosome,  start, end ,community_id, region_length, gc_content, tss_distance, score_phylo))
}


profile_community_biology_hg38 <- function(df_total) {
  
  # 1. Prepare Annotations
  genome <- BSgenome.Hsapiens.UCSC.hg38
  phast <- phastCons100way.UCSC.hg38
  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
  all_tss <- promo_ranges <- promoters(txdb, upstream=0, downstream=1)
  returning_df <- data.frame()
  for ( community_name in unique(df_total$community)){
    df <- df_total |> filter(community == community_name)
    # Convert df to GRanges
    gr <- makeGRangesFromDataFrame(df, keep.extra.columns = TRUE , na.rm = TRUE)
    
    # --- Feature: Region Length ---
    df$region_length <- width(gr)
    
    # --- Feature: GC Content ---
    seqs <- getSeq(genome, gr)
    df$gc_content <- as.numeric(letterFrequency(seqs, letters = "GC", as.prob = TRUE))
    
    # --- Feature: TSS Distance ---
    # Find distance to the nearest Transcription Start Site
    dist_to_tss <- distanceToNearest(gr, all_tss)
    df$tss_distance <- mcols(dist_to_tss)$distance
    
    df$score_phylo <- score(phast, gr)
    # Keep only the new features and the community ID
    df$community_id <- community_name
    returning_df <- rbind(returning_df , df)
  }
  return(returning_df %>% select( chromosome,  start, end ,community_id, region_length, gc_content, tss_distance, score_phylo))
}

load_communties_from_folder <- function(path_comms) {
  df_final <- data.frame()
  files_old_list <- list.files(path = path_comms , full.names = TRUE)
  for(x in seq_along(files_old_list)) {
    df_tmp <- read_delim(files_old_list[[x]] , delim = "\t" , col_names = FALSE) 
    df_tmp$community <- x
    df_final <- rbind(df_final , df_tmp)
  }
  df_final <- df_final |> drop_na()
  colnames(df_final) <- c('chromosome','start','end','community')
  return(df_final)
}
comms_thp1_df <- load_communties_from_folder("/mnt/nas-safu02/sdigiove_workspace/check_th_TCHRONET/new_tchroent_ties/thp1/beds/hg38/")
comms_liver_df <- load_communties_from_folder("/mnt/nas-safu02/sdigiove_workspace/check_th_TCHRONET/new_tchroent_ties/liver/bed/")
comms_ball_df <- load_communties_from_folder("/mnt/nas-safu02/sdigiove_workspace/check_th_TCHRONET/new_tchroent_ties/ball/beds/")


thp1_scores <- profile_community_biology_hg38(comms_thp1_df)
liver_scores <- profile_community_biology_mm10(comms_liver_df)
ball_scores <- profile_community_biology_hg38(comms_ball_df)



##### --- THP1 PLOT AND TEST  ---
thp1_scores$community_id <- as.factor(thp1_scores$community_id)
thp1_scores2 <- thp1_scores |> select(-chromosome , -start, -end)
# Reshape data for plotting
plot_data <- thp1_scores2 %>%
  pivot_longer(cols = -community_id, names_to = "Feature", values_to = "Value")

# 1. Transform only the features that span multiple orders of magnitude
plot_data_transformed <- plot_data %>%
  mutate(Value = case_when(
    Feature %in% c("region_length", "tss_distance") ~ log10(Value + 1),
    TRUE ~ Value
  )) %>%
  # Rename features to reflect the transformation for the labels
  mutate(Feature = case_when(
    Feature == "region_length" ~ "Regions Length (log10)",
    Feature == "tss_distance" ~ "TSS distance (log10)",
    Feature == "gc_contet" ~ "GC Content",
    Feature == "score_phylo" ~ "Phylo Score",
    TRUE ~ Feature
  ))

# 2. Plot as usual
pdf("/mnt/nas-safu02/sdigiove_workspace/check_th_TCHRONET/new_tchroent_ties/thp1/communities_features.pdf"  ,
  height = 5, 
  width = 9,
  family = "ArialMT",
  useDingbats = FALSE)
ggplot(plot_data_transformed, aes(x = community_id, y = Value, fill = community_id)) +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(~Feature, scales = "free_y") +
  theme_classic() +
  theme(
  strip.background = element_blank(), # Removes the box
  panel.border = element_rect(color = "black", fill = NA) # Keeps a border around the data
  )+
  theme (text = element_text(size = 15))+ 
  labs(y = "Transformed Value")
dev.off()


# 1. Perform Pairwise Wilcoxon Tests for all features
stats_df <- plot_data %>%
  group_by(Feature) %>%
  wilcox_test(Value ~ community_id) %>%
  # Add significance stars automatically for the dataframe
  add_significance() %>%
  # Adjust p-values for multiple testing (important for reviewers!)
  adjust_pvalue(method = "fdr")%>%
  add_significance() |> 
  dplyr::select(-.y.)
# 2. Perform Global Test (Kruskal-Wallis)
global_stats <- plot_data %>%
  group_by(Feature) %>%
  kruskal_test(Value ~ community_id)%>%
  add_significance()|> 
  dplyr::select(-.y.)
write.csv(global_stats, "/mnt/nas-safu02/sdigiove_workspace/check_th_TCHRONET/new_tchroent_ties/thp1/thp1_global_statistics_kruskal.csv", quote=F , row.names = FALSE)
write.csv(stats_df, "/mnt/nas-safu02/sdigiove_workspace/check_th_TCHRONET/new_tchroent_ties/thp1/thp1_pairwise_statistics_wilcoxon.csv",quote=F , row.names = FALSE)
write.csv(thp1_scores, "/mnt/nas-safu02/sdigiove_workspace/check_th_TCHRONET/new_tchroent_ties/thp1/Thp1_score_features.csv",quote=F , row.names = FALSE)


##### --- LIVER PLOT AND TEST  ---
liver_scores$community_id <- as.factor(liver_scores$community_id)
liver_scores2 <- liver_scores |> select(-chromosome , -start, -end)
# Reshape data for plotting
plot_data <- liver_scores2 %>%
  pivot_longer(cols = -community_id, names_to = "Feature", values_to = "Value")

# 1. Transform only the features that span multiple orders of magnitude
plot_data_transformed <- plot_data %>%
  mutate(Value = case_when(
    Feature %in% c("region_length", "tss_distance") ~ log10(Value + 1),
    TRUE ~ Value
  )) %>%
  # Rename features to reflect the transformation for the labels
  mutate(Feature = case_when(
    Feature == "region_length" ~ "Regions Length (log10)",
    Feature == "tss_distance" ~ "TSS distance (log10)",
    Feature == "gc_contet" ~ "GC Content",
    Feature == "score_phylo" ~ "Phylo Score",
    TRUE ~ Feature
  ))

# 2. Plot as usual
pdf("/mnt/nas-safu02/sdigiove_workspace/check_th_TCHRONET/new_tchroent_ties/liver/communities_features.pdf"  ,
  height = 5, 
  width = 9,
  family = "ArialMT",
  useDingbats = FALSE)
ggplot(plot_data_transformed, aes(x = community_id, y = Value, fill = community_id)) +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(~Feature, scales = "free_y") +
  theme_classic() +
  theme(
  strip.background = element_blank(), # Removes the box
  panel.border = element_rect(color = "black", fill = NA) # Keeps a border around the data
  )+
  theme (text = element_text(size = 15))+ 
  labs(y = "Transformed Value")
dev.off()


# 1. Perform Pairwise Wilcoxon Tests for all features
stats_df <- plot_data %>%
  group_by(Feature) %>%
  wilcox_test(Value ~ community_id) %>%
  # Add significance stars automatically for the dataframe
  add_significance() %>%
  # Adjust p-values for multiple testing (important for reviewers!)
  adjust_pvalue(method = "fdr")%>%
  add_significance() |> 
  dplyr::select(-.y.)
# 2. Perform Global Test (Kruskal-Wallis)
global_stats <- plot_data %>%
  group_by(Feature) %>%
  kruskal_test(Value ~ community_id)%>%
  add_significance()|> 
  dplyr::select(-.y.)
write.csv(global_stats, "/mnt/nas-safu02/sdigiove_workspace/check_th_TCHRONET/new_tchroent_ties/liver/liver_global_statistics_kruskal.csv", quote=F , row.names = FALSE)
write.csv(stats_df, "/mnt/nas-safu02/sdigiove_workspace/check_th_TCHRONET/new_tchroent_ties/liver/liver_pairwise_statistics_wilcoxon.csv",quote=F , row.names = FALSE)
write.csv(liver_scores, "/mnt/nas-safu02/sdigiove_workspace/check_th_TCHRONET/new_tchroent_ties/liver/Liver_score_features.csv",quote=F , row.names = FALSE)



##### --- BALL PLOT AND TEST  ---
ball_scores$community_id <- as.factor(ball_scores$community_id)
ball_scores2 <- ball_scores |> select(-chromosome , -start, -end)
# Reshape data for plotting
plot_data <- ball_scores2 %>%
  pivot_longer(cols = -community_id, names_to = "Feature", values_to = "Value")

# 1. Transform only the features that span multiple orders of magnitude
plot_data_transformed <- plot_data %>%
  mutate(Value = case_when(
    Feature %in% c("region_length", "tss_distance") ~ log10(Value + 1),
    TRUE ~ Value
  )) %>%
  # Rename features to reflect the transformation for the labels
  mutate(Feature = case_when(
    Feature == "region_length" ~ "Regions Length (log10)",
    Feature == "tss_distance" ~ "TSS distance (log10)",
    Feature == "gc_contet" ~ "GC Content",
    Feature == "score_phylo" ~ "Phylo Score",
    TRUE ~ Feature
  ))

# 2. Plot as usual
pdf("/mnt/nas-safu02/sdigiove_workspace/check_th_TCHRONET/new_tchroent_ties/ball/communities_features.pdf"  ,
  height = 5, 
  width = 9,
  family = "ArialMT",
  useDingbats = FALSE)
ggplot(plot_data_transformed, aes(x = community_id, y = Value, fill = community_id)) +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(~Feature, scales = "free_y") +
  theme_classic() +
  theme(
  strip.background = element_blank(), # Removes the box
  panel.border = element_rect(color = "black", fill = NA) # Keeps a border around the data
  )+
  theme (text = element_text(size = 15))+ 
  labs(y = "Transformed Value")
dev.off()


# 1. Perform Pairwise Wilcoxon Tests for all features
stats_df <- plot_data %>%
  group_by(Feature) %>%
  wilcox_test(Value ~ community_id) %>%
  # Add significance stars automatically for the dataframe
  add_significance() %>%
  # Adjust p-values for multiple testing (important for reviewers!)
  adjust_pvalue(method = "fdr")%>%
  add_significance() |> 
  dplyr::select(-.y.)
# 2. Perform Global Test (Kruskal-Wallis)
global_stats <- plot_data %>%
  group_by(Feature) %>%
  kruskal_test(Value ~ community_id)%>%
  add_significance()|> 
  dplyr::select(-.y.)
write.csv(global_stats, "/mnt/nas-safu02/sdigiove_workspace/check_th_TCHRONET/new_tchroent_ties/ball/ball_global_statistics_kruskal.csv", quote=F , row.names = FALSE)
write.csv(stats_df, "/mnt/nas-safu02/sdigiove_workspace/check_th_TCHRONET/new_tchroent_ties/ball/ball_pairwise_statistics_wilcoxon.csv",quote=F , row.names = FALSE)
write.csv(ball_scores, "/mnt/nas-safu02/sdigiove_workspace/check_th_TCHRONET/new_tchroent_ties/ball/Ball_score_features.csv",quote=F , row.names = FALSE)
