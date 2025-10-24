suppressMessages (expr = { 
    library(rGREAT)
    library(tidyverse)
    library(clustree)
    library(GSEABase)
    library(clusterProfiler)
    library(enrichR)
    library(msigdbr)
})

cluster_df  <- read_delim("~/T-ChroNet/paper_analysis/data/THP1/results/cluster_tree_all_res.tsv" , delim = "\t")
png("~/T-ChroNet/paper_analysis/data/THP1/results/cluster_tree.png" , height = 5 , width = 7 , units = 'in' , res = 300)
clustree(cluster_df, prefix = "Cluster_")
dev.off()

cluster_df  <- read_delim("~/T-ChroNet/paper_analysis/data/LiverDevelopment/results/cluster_tree_all_res.tsv" , delim = "\t") |>
select(-'Cluster_1.81',-'Cluster_1.905',-'Cluster_2.0')
#png("~/T-ChroNet/paper_analysis/data/LiverDevelopment/results/clustertree.png" , height = 5 , width = 7 , units = 'in' , res = 300)
clustree(cluster_df, prefix = "Cluster_")
#dev.off()

cluster_df  <- read_delim("~/T-ChroNet/paper_analysis/data/BCP-ALL/results/cluster_tree_all_res_th05_nfcore.tsv" , delim = "\t")
#png("~/T-ChroNet/paper_analysis/data/BCP-ALL/results/clustertree.png" , height = 5 , width = 7 , units = 'in' , res = 300)
clustree(cluster_df, prefix = "Cluster_")
#dev.off()