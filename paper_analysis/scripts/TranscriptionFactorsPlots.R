library(tidyverse)
library(ComplexHeatmap)

# THP-1
folders_to_compare = list.files('/mnt/nas-safu/analysis/PhDsdigiove/method_coAcces/data/CellReport/communities_k27/res_1_5/cistrom/')

for (community in folders_to_compare){ 
    to_read = paste('/mnt/nas-safu/analysis/PhDsdigiove/method_coAcces/data/CellReport/communities_k27/res_1_5/cistrom/' ,community , sep = "")

    homer_tfs = read_delim(to_read , delim = ",")

    #homer_tfs  <- homer_tfs |> filter(grepl('liver' , Biosource ) | grepl('Liver' , Biosource ))

    homer_tfs = homer_tfs |> group_by(Factor) |> summarise(max_value = max(GIGGLE_score)) |> select(Factor , max_value)
    

    colnames(homer_tfs) = c('Factor' , community)

    if (community == folders_to_compare[1]){
        final_homer_tfs = homer_tfs}
    else {
        final_homer_tfs = merge.data.frame(final_homer_tfs , homer_tfs , all = T)}
}
final_homer_tfs <- final_homer_tfs |> column_to_rownames("Factor")

tfs_name = read_delim("/mnt/nas-safu/analysis/PhDsdigiove/method_coAcces/data/TF_names_v_1.01.txt" , delim = "\t" , col_names = F)
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
theme(legend.position="top" , legend.text = element_text(size=10) ) +
  #axis.text.y = element_text(color = "grey20", size = 10, angle = 0, hjust = 1, vjust = 0) ) +
scale_x_discrete( limits = unique(final_homer_scatterplot$Community) )
# theme(axis.title.y=element_blank(), #axis.text.y = element_text(color = "grey20", size = 10, angle = 0, hjust = 1, vjust = 0, face = "plain") ,
#       axis.text.y=element_blank(),
#       legend.position="none")
ggsave("/mnt/nas-safu/analysis/PhDsdigiove/method_coAcces/data/CellReport/pictures/TFs_R_withlabels_THP1.png", width = 7, height = 8, units = "in", dpi = 300)
# ggsave("/mnt/nas-safu/analysis/PhDsdigiove/method_coAcces/data/CellReport/pictures/TFs_R_Nolabels_THP1.png", width = 7, height = 8, units = "in", dpi = 300)


# LIVER DEVELOPMENT
folders_to_compare = list.files('/mnt/nas-safu/analysis/PhDsdigiove/method_coAcces/data/liver_devel_mouse_ENCODE/communities/res_0_8/cistromdb/')

for (community in folders_to_compare){ 
    to_read = paste('/mnt/nas-safu/analysis/PhDsdigiove/method_coAcces/data/liver_devel_mouse_ENCODE/communities/res_0_8/cistromdb/' ,community , sep = "")

    homer_tfs = read_delim(to_read , delim = ",")

    #homer_tfs  <- homer_tfs |> filter(grepl('liver' , Biosource ) | grepl('Liver' , Biosource ))

    homer_tfs = homer_tfs |> group_by(Factor) |> summarise(max_value = max(GIGGLE_score)) |> select(Factor , max_value)
    

    colnames(homer_tfs) = c('Factor' , community)

    if (community == folders_to_compare[1]){
        final_homer_tfs = homer_tfs}
    else {
        final_homer_tfs = merge.data.frame(final_homer_tfs , homer_tfs , all = T)}
}
final_homer_tfs <- final_homer_tfs |> column_to_rownames("Factor")

tfs_name = read_delim("/mnt/nas-safu/analysis/PhDsdigiove/method_coAcces/data/mouse_tfs.txt" , delim = "\t")
tfs_name <- tfs_name %>%
  mutate(across(where(is.character), toupper))
real_tfs  <- c()
for (x in rownames(final_homer_tfs)){
    if (x %in% tfs_name[['Tfs']]){
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
  theme(legend.position="top" , legend.text = element_text(size=10) ,
  axis.text.y = element_text(color = "grey20", size = 10, angle = 0, hjust = 1, vjust = 0) ) 
  #scale_x_discrete(limits =colnames(clusteing_df), guide = guide_axis(angle = 90)) +
  # theme(axis.title.y=element_blank(), #axis.text.y = element_text(color = "grey20", size = 10, angle = 0, hjust = 1, vjust = 0, face = "plain") ,
  #       axis.text.y=element_blank(),
  #       legend.position="none")
# ggsave("/mnt/nas-safu/analysis/PhDsdigiove/method_coAcces/data/liver_devel_mouse_ENCODE/communities/res_0_8/TFs_R_withlabels_liver.png", width = 7, height = 8, units = "in", dpi = 300)
# ggsave("/mnt/nas-safu/analysis/PhDsdigiove/method_coAcces/data/liver_devel_mouse_ENCODE/communities/res_0_8/TFs_R_Nolabels_liver.png", width = 7, height = 5, units = "in", dpi = 300)



# B-ALL
folders_to_compare = list.files('/mnt/nas-safu/analysis/PhDsdigiove/method_coAcces/data/BALL/communities_broad/res_0_9/cistromdf_broad/')

for (community in folders_to_compare){ 
    to_read = paste('/mnt/nas-safu/analysis/PhDsdigiove/method_coAcces/data/BALL/communities_broad/res_0_9/cistromdf_broad/' ,community , sep = "")

    homer_tfs = read_delim(to_read , delim = ",")

    #homer_tfs  <- homer_tfs |> filter(grepl('liver' , Biosource ) | grepl('Liver' , Biosource ))

    homer_tfs = homer_tfs |> group_by(Factor) |> summarise(max_value = max(GIGGLE_score)) |> select(Factor , max_value)
    

    colnames(homer_tfs) = c('Factor' , community)

    if (community == folders_to_compare[1]){
        final_homer_tfs = homer_tfs}
    else {
        final_homer_tfs = merge.data.frame(final_homer_tfs , homer_tfs , all = T)}
}
final_homer_tfs <- final_homer_tfs |> column_to_rownames("Factor")

tfs_name = read_delim("/mnt/nas-safu/analysis/PhDsdigiove/method_coAcces/data/TF_names_v_1.01.txt" , delim = "\t" , col_names = F)
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
final_homer_scatterplot$Community <- as.numeric(sub(".*?(\\d).*", "\\1", final_homer_scatterplot$Community))

ggplot(data = final_homer_scatterplot , aes(x = Community, y = Factor ) ) +
  geom_point( aes( size = ComboScore) ,  shape = 21  , stroke = 0.6  , color = 'black') +
  geom_point( aes( size = ComboScore) ,  shape = 21  , fill = "black" , stroke = 0.6 , alpha= 0.4 , color = 'black') +
  theme_classic() +
  ylab("")+
    theme(legend.position="top" , legend.text = element_text(size=10) ,
    axis.text.y = element_text(color = "grey20", size = 10, angle = 0, hjust = 1, vjust = 0) )  +
    scale_x_discrete( limits = unique(final_homer_scatterplot$Community) )
    # theme(axis.title.y=element_blank(), #axis.text.y = element_text(color = "grey20", size = 10, angle = 0, hjust = 1, vjust = 0, face = "plain") ,
    #       axis.text.y=element_blank(),
    #       legend.position="none")
# ggsave("/mnt/nas-safu/analysis/PhDsdigiove/method_coAcces/data/BALL/communities_broad/res_0_9/TFs_R_withlabels_BCPALL.png", width = 7, height = 8, units = "in", dpi = 300)
# ggsave("/mnt/nas-safu/analysis/PhDsdigiove/method_coAcces/data/liver_devel_mouse_ENCODE/communities/res_0_8/TFs_R_Nolabels_liver.png", width = 7, height = 5, units = "in", dpi = 300)