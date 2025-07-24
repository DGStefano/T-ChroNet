suppressMessages (expr = { 
    library(rGREAT)
    library(tidyverse)
    library(clustree)
    library(GSEABase)
    library(clusterProfiler)
    library(enrichR)
    library(msigdbr)
})

cluster_df  <- read_delim("/mnt/nas-safu/analysis/PhDsdigiove/method_coAcces/data/CellReport/communities_k27/res_1_5/cluster_tree_all_res.tsv" , delim = "\t")
png("/mnt/nas-safu/analysis/PhDsdigiove/method_coAcces/pictures/cell_report/clustertree.png" , height = 5 , width = 7 , units = 'in' , res = 300)
clustree(cluster_df, prefix = "Cluster_")
dev.off()

cluster_df  <- read_delim("/mnt/nas-safu/analysis/PhDsdigiove/method_coAcces/data/liver_devel_mouse_ENCODE/communities/cluster_tree_all_res.tsv" , delim = "\t") |>
select(-'Cluster_1.81',-'Cluster_1.905',-'Cluster_2.0')
#png("/mnt/nas-safu/analysis/PhDsdigiove/method_coAcces/pictures/cell_report/clustertree.png" , height = 5 , width = 7 , units = 'in' , res = 300)
clustree(cluster_df, prefix = "Cluster_")
#dev.off()

cluster_df  <- read_delim("/mnt/nas-safu/analysis/PhDsdigiove/method_coAcces/data/liver_devel_mouse_ENCODE/communities/cluster_tree_all_res_th05_nfcore.tsv" , delim = "\t")
#png("/mnt/nas-safu/analysis/PhDsdigiove/method_coAcces/pictures/cell_report/clustertree.png" , height = 5 , width = 7 , units = 'in' , res = 300)
clustree(cluster_df, prefix = "Cluster_")
#dev.off()