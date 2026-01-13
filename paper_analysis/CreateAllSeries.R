library(clustree)
library(tidyverse)
library(readr)
library(tidyverse)
library(ComplexHeatmap)
library(clusterProfiler)
library(msigdbr)
library(TChroNetR)
# edge_files <- list.files("/mnt/nas-safu02/sdigiove_workspace/check_th_TCHRONET/new_tchroent_ties/thp1/th/" , full.names = T)
# matrix_path <- "~/T-ChroNet/paper_analysis/data/THP1/lognorm_edgeR_limma_countsInCellReport_TheiCounts_NoStatic_mean.tsv"
# series_thp1 <- build_TChroNetSeries_object(edge_files,matrix_path ,method = "Leiden",resolutions_list = seq(0.5, 1.5, 0.1),run_cd = T,seed = 1234 , transitivity = FALSE)
# write_rds(series_thp1 , "/mnt/nas-safu02/sdigiove_workspace/check_th_TCHRONET/new_tchroent_ties/thp1/series_obj_thp1.rds")

# edge_files_liver <- list.files("/mnt/nas-safu02/sdigiove_workspace/check_th_TCHRONET/new_tchroent_ties/liver/th/" , full.names = T)
# matrix_path <- "~/T-ChroNet/paper_analysis/data/LiverDevelopment/normalized_samplemean_multicov_all_sites_all_timepoint.tsv"
# seiries_obj_liver <- build_TChroNetSeries_object(edge_files_liver,matrix_path,method = "Leiden",resolutions_list = seq(0.5, 1.5, 0.1),run_cd = T,seed = 123, transitivity = FALSE)
# write_rds(seiries_obj_liver , "/mnt/nas-safu02/sdigiove_workspace/check_th_TCHRONET/new_tchroent_ties/liver/series_obj_liver.rds")

edge_files_ball <- list.files("/mnt/nas-safu02/sdigiove_workspace/check_th_TCHRONET/new_tchroent_ties/ball/th/" , full.names = T)
matrix_path <- "~/T-ChroNet/paper_analysis/data/BCP-ALL/allPatientsPatientsNOMERGED_Ball_Multicov_nfcore.tsv"
series_ball <- build_TChroNetSeries_object(edge_files_ball,matrix_path,method = "Leiden",resolutions_list = seq(0.5, 1.5, 0.1),run_cd = T,seed = 123, transitivity = FALSE)
write_rds(series_ball , "/mnt/nas-safu02/sdigiove_workspace/check_th_TCHRONET/new_tchroent_ties/ball/series_obj_ball.rds")
