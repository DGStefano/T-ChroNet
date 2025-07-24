library(tidyverse)
library(edgeR)
library(Matrix)
library(readxl)
library(pheatmap)
library(limma)
library(edgeR)
library(MatrixGenerics)
library(tidyverse)


# FUNCTIONS
normlization_cpm = function (df ,groups) {
  e=data.matrix(df)
  e=subset(e, rowMeans(e) > 5)
  e=subset(e, rowMeans(e) < 5000)
  
  #make DGE list object(EdgeR basic object)
  regions <- DGEList(counts=e, group=groups)
  
  #check if opportune correct for the common dispersion by using housekeeping genes as suggested on EDGEr manual (do not do)
  keep <- filterByExpr(regions , group=groups , min.prop = 0.1)
  regions<- regions[keep,,keep.lib.sizes=FALSE]
  
  #EdgeR data correction
  regions <- estimateCommonDisp(regions)
  regions <- estimateTagwiseDisp(regions)
  
  #TMM normalization
  regions <- calcNormFactors(regions, method="TMM" )
  cpm_tmm <- cpm(regions , log = TRUE) #
  return( list("regions" = regions , "cpm" = cpm_tmm ))
}

# Function to select genes present in both replicates for each sample
filter_genes <- function(data) {
  # Find replicate pairs by sample (assumes naming convention "SampleX_RepY")
  replicate_pairs <- split(colnames(data), sub("_REP\\d+$", "", colnames(data)))

  # Logical vector for genes present in all replicate pairs
  # keep_genes <- Reduce(`&`, lapply(replicate_pairs, function(replicates) {
  #   # Ensure the genes are expressed in all replicates for the sample
  #   rowSums(data[, replicates]) == length(replicates)
  # }))


  # Logical vector for peaks present in both replicates for at least one condition
  keep_peaks <- Reduce(`|`, lapply(replicate_pairs, function(replicates) {
    # Check if the peak is present (value != 0) in both replicates for the sample
    rowSums(data[, replicates] != -3.474) == length(replicates)
  }))
  # Filter the data for these genes
  data[keep_peaks, ]
}



# THP1 NORMALIZATION

counts <- read.delim("/mnt/nas-safu/analysis/PhDsdigiove/method_coAcces/data/CellReport/out_mutlicov_matrix/Multicov_AllAnalyzedPeaksWithAnnotations_NOSTATIC.tsv", sep = "\t") |> 
  unite('peaks', c('chromosome', 'start', 'end'), sep = "-")
counts <- counts[!duplicated(counts$peaks),] |> rownames_to_column("_") |> select(-'_') |> column_to_rownames('peaks')

counts <- read.delim("/mnt/nas-safu/analysis/PhDsdigiove/method_coAcces/data/CellReport/their_regions_counts.bed" , #
                     sep = "\t" ) |> unite('peaks' , c('chr','start','end') , sep ="-")
counts <- counts[!duplicated(counts$peaks),] |> rownames_to_column("_") |> select(-'_') |> column_to_rownames('peaks')
sample.annot = c('t0','t0','t30min','t30min','t60min','t60min','t90min','t90min','t120min','t120min','t240min' ,'t240min' , 't360min', 't360min', 't1440min', 't1440min') 

replicates <- rep(c('1','2') , times = 8)

# Make a DGEList and add metadata
y <- DGEList(counts = counts)

y$samples$timepoint <- factor(sample.annot)
y$samples$replicated <- factor(replicates)

# Normalize and obtain logcounts for QC
y <- calcNormFactors(y)
logCPMs <- cpm(y, log = TRUE)


# Run PCA
pca <- prcomp(t(logCPMs))

# Combine PCA coordinates with the metadata from the DGEList
to_plot <- data.frame(pca$x, y$samples)

# Calculate how many % of total variance is explained by each principal component
percentVar <- pca$sdev^2 / sum( pca$sdev^2 )*100

# We focus here on PC1 and PC2
use.pcs <- c(1,2)
labs <- paste0(paste0("PC", use.pcs, " - "), paste0("Var.expl = ", round(percentVar[use.pcs], 2), "%"))

ggplot(to_plot, aes(x=PC1, y=PC2, shape=replicates, color=timepoint)) + 
  geom_point(size=3) +
  xlab(labs[1]) + ylab(labs[2])
# ggsave("/mnt/nas-safu/analysis/PhDsdigiove/method_coAcces/pictures/cell_report/batch_correction/BEFORE_correction_PCA.png" , height = 5 , width = 9 , units = 'in' , dpi = 300)

# correct for the batch which here is the "kit"
batch <- factor(y$samples$replicated)
timepoints <- factor(y$samples$timepoint)

logCPMs_corrected <- limma::removeBatchEffect(logCPMs, batch = batch , group = timepoints)

# Run PCA
pca <- prcomp(t(logCPMs_corrected))

# Combine PCA coordinates with the metadata from the DGEList
to_plot <- data.frame(pca$x, y$samples)

# Calculate how many % of total variance is explained by each principal component
percentVar <- pca$sdev^2 / sum( pca$sdev^2 )*100

# We focus here on PC1 and PC2
use.pcs <- c(1,2)
labs <- paste0(paste0("PC", use.pcs, " - "), paste0("Var.expl = ", round(percentVar[use.pcs], 2), "%"))

ggplot(to_plot, aes(x=PC1, y=PC2, shape=replicates, color=timepoint)) + 
  geom_point(size=3) +
  xlab(labs[1]) + ylab(labs[2])
# ggsave("/mnt/nas-safu/analysis/PhDsdigiove/method_coAcces/pictures/cell_report/batch_correction/AFTER_correction_PCA.png" , height = 5 , width = 9 , units = 'in' , dpi = 300)

logCPMs_corrected_filtered <- logCPMs_corrected#[rownames(counts),]

cor_matrix <- cor(logCPMs_corrected_filtered , method = "pearson") # Use 'spearman' or 'kendall' if needed


# Plot the heatmap
# png("/mnt/nas-safu/analysis/PhDsdigiove/method_coAcces/pictures/cell_report/batch_correction/AFTER_correction_PCA.png" , height = 5 , width = 9 , units = 'in' , res = 300)
pheatmap(cor_matrix, 
         display_numbers = TRUE, # Show correlation values on the heatmap
         clustering_method = "complete", # Clustering method
         color = colorRampPalette(c("blue", "white", "red"))(100)) # Color gradient
# dev.off()

row_mean = rowMeans2(logCPMs_corrected_filtered)
row_var = rowVars(logCPMs_corrected_filtered)

data.frame("rowMean" = row_mean , "rowVar" = row_var ) |> 
  ggplot(aes(x = rowMean , y = rowVar))+
  geom_point() +
  geom_hline(yintercept = 1.0)
# ggsave("/mnt/nas-safu/analysis/PhDsdigiove/method_coAcces/pictures/cell_report/batch_correction/Means_vs_Variance_TH1.png" , height = 5 , width = 9 , units = 'in' , dpi = 300)

dim(logCPMs_corrected_filtered)

# as.data.frame(logCPMs_corrected_filtered) |> rownames_to_column('peaks') |> write_delim("/mnt/nas-safu/analysis/PhDsdigiove/method_coAcces/data/CellReport/out_mutlicov_matrix/lognorm_edgeR_limma_countsInCellReport_TheiCounts_All.tsv" , delim = "\t" )

# LIVER DEVELOPMENT
#### STARDARD EDGER AND LIMMA ANALYSIS
# counts <- read.delim("/mnt/nas-safu/analysis/PhDsdigiove/method_coAcces/data/heartdevel/multicov_mouse_accessibility_allasterlist.txt" , #
#                      sep = "\t" ) |> unite('peaks' , c('chromosome','start','end') , sep ="-")
counts <- read.delim("/mnt/nas-safu/analysis/PhDsdigiove/method_coAcces/data/liver_devel_mouse_ENCODE/multicov_all_sites_all_timepoint.tsv" , #
                     sep = "\t" ) |> unite('peaks' , c('chromosome','start','end') , sep ="-")                     
counts <- counts[!duplicated(counts$peaks),] |> rownames_to_column("_") |> dplyr::select(-'_') |> column_to_rownames('peaks')


#removing rep 2 tp 13 due to high discrepacy with rep 1
# counts <- counts |> dplyr::select(-musmusculus_heart_13_5d_rep2 , -musmusculus_heart_13_5d_rep1)

sample.annot = c('t11_5','t11_5' , 't12_5','t12_5','t13_5','t13_5' , 't14_5','t14_5','t15_5','t15_5' , 't16_5','t16_5')  #

replicates <- rep(c('1','2') , times = 6)
# replicates <- rep('1','2','1','2','1','1','2','1','2','1','2' )

# Make a DGEList and add metadata
y <- DGEList(counts = counts)

y$samples$timepoint <- factor(sample.annot)
y$samples$replicated <- factor(replicates)

# Normalize and obtain logcounts for QC
y <- calcNormFactors(y)
logCPMs <- cpm(y, log = TRUE)


# Run PCA
pca <- prcomp(t(logCPMs))

# Combine PCA coordinates with the metadata from the DGEList
to_plot <- data.frame(pca$x, y$samples)

# Calculate how many % of total variance is explained by each principal component
percentVar <- pca$sdev^2 / sum( pca$sdev^2 )*100

# We focus here on PC1 and PC2
use.pcs <- c(1,2)
labs <- paste0(paste0("PC", use.pcs, " - "), paste0("Var.expl = ", round(percentVar[use.pcs], 2), "%"))

ggplot(to_plot, aes(x=PC1, y=PC2, shape=replicates, color=timepoint)) + 
  geom_point(size=3) +
  xlab(labs[1]) + ylab(labs[2])
ggsave("/mnt/nas-safu/analysis/PhDsdigiove/method_coAcces/data/liver_devel_mouse_ENCODE/pictures/BEFORE_correction_PCA.png" , height = 5 , width = 9 , units = 'in' , dpi = 300)

# correct for the batch which here is the "kit"
batch <- factor(y$samples$replicated)
timepoints <- factor(y$samples$timepoint)

logCPMs_corrected <- limma::removeBatchEffect(logCPMs, group = timepoints)
# logCPMs_corrected  <- logCPMs

# Run PCA
pca <- prcomp(t(logCPMs_corrected))

# Combine PCA coordinates with the metadata from the DGEList
to_plot <- data.frame(pca$x, y$samples)

# Calculate how many % of total variance is explained by each principal component
percentVar <- pca$sdev^2 / sum( pca$sdev^2 )*100

# We focus here on PC1 and PC2
use.pcs <- c(1,2)
labs <- paste0(paste0("PC", use.pcs, " - "), paste0("Var.expl = ", round(percentVar[use.pcs], 2), "%"))

ggplot(to_plot, aes(x=PC1, y=PC2, shape=replicates, color=timepoint)) + 
  geom_point(size=3) +
  xlab(labs[1]) + ylab(labs[2])
ggsave("/mnt/nas-safu/analysis/PhDsdigiove/method_coAcces/data/liver_devel_mouse_ENCODE/pictures/AFTER_correction_PCA.png" , height = 5 , width = 9 , units = 'in' , dpi = 300)

logCPMs_corrected_filtered <- logCPMs_corrected #|> as.data.frame()

#logCPMs_corrected_filtered$musmusculus_heart_13_5d_rep2  <- logCPMs_corrected_filtered$musmusculus_heart_13_5d_rep1

cor_matrix <- cor(logCPMs_corrected_filtered , method = "pearson") # Use 'spearman' or 'kendall' if needed


# Plot the heatmap
ph <- pheatmap(cor_matrix, 
         display_numbers = TRUE, # Show correlation values on the heatmap
         clustering_method = "complete", # Clustering method
         color = colorRampPalette(c("blue", "white", "red"))(100)) # Color gradient
#png("/mnt/nas-safu/analysis/PhDsdigiove/method_coAcces/data/liver_devel_mouse_ENCODE/pictures/AFTER_correction_Correlation.png" , height = 5 , width = 9 , units = 'in' , res = 300)
ph
#dev.off()

logCPMs_corrected_filtered_noextra <- logCPMs_corrected_filtered[!grepl("chrX", rownames(logCPMs_corrected_filtered)) , ]
logCPMs_corrected_filtered_noextra <- logCPMs_corrected_filtered_noextra[!grepl("chrY", rownames(logCPMs_corrected_filtered_noextra)) , ]
logCPMs_corrected_filtered_noextra <- logCPMs_corrected_filtered_noextra[!grepl("gl0", rownames(logCPMs_corrected_filtered_noextra)) , ]
logCPMs_corrected_filtered_noextra <- logCPMs_corrected_filtered_noextra[!grepl("random", rownames(logCPMs_corrected_filtered_noextra)) , ]
logCPMs_corrected_filtered_noextra <- logCPMs_corrected_filtered_noextra[!grepl("chrM", rownames(logCPMs_corrected_filtered_noextra)) , ]

row_mean = rowMeans2(logCPMs_corrected_filtered_noextra)
row_var = rowVars(logCPMs_corrected_filtered_noextra)

data.frame("rowMean" = row_mean , "rowVar" = row_var ) |> 
  ggplot(aes(x = rowMean , y = rowVar))+
  geom_point() +
  geom_hline(yintercept = 0.1)


logCPMs_var  <- logCPMs_corrected_filtered_noextra[row_var > 0.1,] |> as.data.frame()

result <- data.frame(
  lapply(seq(1, ncol(logCPMs_var), by = 2), function(i) {
    rowMeans(logCPMs_var[, i:min(i + 1, ncol(logCPMs_var))])
  })
)

colnames(result)  <- c('t11_5', 't12_5','t13_5' , 't14_5','t15_5', 't16_5') 

# as.data.frame(result) |> rownames_to_column("peaks") |> write_delim("/mnt/nas-safu/analysis/PhDsdigiove/method_coAcces/data/liver_devel_mouse_ENCODE/network_files/normalized_samplemean_multicov_all_sites_all_timepoint.tsv" , delim = "\t")

# B-ALL
#### STARDARD EDGER AND LIMMA ANALYSIS
# counts <- read.delim("/mnt/nas-safu/analysis/PhDsdigiove/method_coAcces/data/BALL/multicovallsampe.bed" , #
#                      sep = "\t" )

counts <- read_delim("/mnt/nas-safu/analysis/PhDsdigiove/method_coAcces/data/BALL/preliminary_analysis/results_atac/bwa/merged_library/macs2/broad_peak/consensus/consensus_peaks.mLb.clN.featureCounts.txt" , #
                     delim = "\t" , comment = "#") |> select(-Geneid , -Strand , -Length)  |> unite('peaks' , c('Chr','Start','End') , sep ="-")
counts <- counts[!duplicated(counts$peaks),] |> rownames_to_column("_") |> dplyr::select(-'_') |> column_to_rownames('peaks')
counts <- counts |> select('patient_6_Healthy_ATACseq_REP1.mLb.clN.sorted.bam','patient_9_Healthy_ATACseq_REP1.mLb.clN.sorted.bam','patient_4_Healthy_ATACseq_REP1.mLb.clN.sorted.bam','patient_7_Healthy_ATACseq_REP1.mLb.clN.sorted.bam','patient_5_Healthy_ATACseq_REP1.mLb.clN.sorted.bam','patient_8_Healthy_ATACseq_REP1.mLb.clN.sorted.bam',
'patient_14_Primary_ATACseq_REP1.mLb.clN.sorted.bam','patient_12_Primary_ATACseq_REP1.mLb.clN.sorted.bam','patient_10_Primary_ATACseq_REP1.mLb.clN.sorted.bam','patient_13_Primary_ATACseq_REP1.mLb.clN.sorted.bam','patient_11_Primary_ATACseq_REP1.mLb.clN.sorted.bam','patient_16_Primary_ATACseq_REP1.mLb.clN.sorted.bam','patient_7_Primary_ATACseq_REP1.mLb.clN.sorted.bam',,'patient_19_Primary_ATACseq_REP1.mLb.clN.sorted.bam','patient_17_Primary_ATACseq_REP1.mLb.clN.sorted.bam','patient_27_Primary_ATACseq_REP1.mLb.clN.sorted.bam','patient_26_Primary_ATACseq_REP1.mLb.clN.sorted.bam',
'patient_21_Remission_ATACseq_REP1.mLb.clN.sorted.bam','patient_22_Remission_ATACseq_REP1.mLb.clN.sorted.bam','patient_18_Remission_ATACseq_REP1.mLb.clN.sorted.bam','patient_26_Remission_ATACseq_REP1.mLb.clN.sorted.bam','patient_20_Remission_ATACseq_REP1.mLb.clN.sorted.bam','patient_25_Remission_ATACseq_REP1.mLb.clN.sorted.bam','patient_23_Remission_ATACseq_REP1.mLb.clN.sorted.bam','patient_24_Remission_ATACseq_REP1.mLb.clN.sorted.bam',
'patient_31_Relapse_ATACseq_REP1.mLb.clN.sorted.bam','patient_28_Relapse_ATACseq_REP1.mLb.clN.sorted.bam','patient_29_Relapse_ATACseq_REP1.mLb.clN.sorted.bam','patient_33_Relapse_ATACseq_REP1.mLb.clN.sorted.bam','patient_32_Relapse_ATACseq_REP1.mLb.clN.sorted.bam','patient_27_Relapse_ATACseq_REP1.mLb.clN.sorted.bam','patient_30_Relapse_ATACseq_REP1.mLb.clN.sorted.bam','patient_14_Relapse_ATACseq_REP1.mLb.clN.sorted.bam')

colnames(counts)  <- c('patient_6_Healthy','patient_9_Healthy','patient_4_Healthy','patient_7_Healthy','patient_5_Healthy','patient_8_Healthy','patient_14_Primary','patient_12_Primary','patient_10_Primary','patient_13_Primary','patient_11_Primary','patient_16_Primary','patient_7_Primary','patient_19_Primary','patient_17_Primary','patient_27_Primary','patient_26_Primary','patient_21_Remission','patient_22_Remission','patient_18_Remission','patient_26_Remission','patient_20_Remission','patient_25_Remission','patient_23_Remission','patient_24_Remission','patient_31_Relapse','patient_28_Relapse','patient_29_Relapse','patient_33_Relapse','patient_32_Relapse','patient_27_Relapse','patient_30_Relapse','patient_14_Relapse')

#removing rep 2 tp 13 due to high discrepacy with rep 1
# counts <- counts |> dplyr::select(-musmusculus_heart_13_5d_rep2 , -musmusculus_heart_13_5d_rep1)

sample.annot = rep(c("healthy" ,"Primary" , "Remission" , "Relapse" ) , c(6, 11 , 8 , 8) ) 


#replicates <- rep(c('1','2') , times = 6)
# replicates <- rep('1','2','1','2','1','1','2','1','2','1','2' )

# Make a DGEList and add metadata
y <- DGEList(counts = counts , group = sample.annot)
keep <- filterByExpr(y)
y <- y[keep,,keep.lib.sizes=FALSE]
y <- normLibSizes(y)
y$samples$timepoint <- factor(sample.annot)
#y$samples$replicated <- factor(replicates)

# Normalize and obtain logcounts for QC
y <- calcNormFactors(y)
logCPMs <- cpm(y, log = TRUE)


# Run PCA
pca <- prcomp(t(logCPMs))

# Combine PCA coordinates with the metadata from the DGEList
to_plot <- data.frame(pca$x, y$samples)

# Calculate how many % of total variance is explained by each principal component
percentVar <- pca$sdev^2 / sum( pca$sdev^2 )*100

# We focus here on PC1 and PC2
use.pcs <- c(1,2)
labs <- paste0(paste0("PC", use.pcs, " - "), paste0("Var.expl = ", round(percentVar[use.pcs], 2), "%"))

ggplot(to_plot, aes(x=PC1, y=PC2,  color=timepoint)) + #,
  geom_point(size=3) +
  xlab(labs[1]) + ylab(labs[2])
# ggsave("/mnt/nas-safu/analysis/PhDsdigiove/method_coAcces/data/BALL/communities_broad/BEFORE_correction_PCA.png" , height = 5 , width = 9 , units = 'in' , dpi = 300)

# correct for the batch which here is the "kit"
#batch <- factor(y$samples$replicated)
timepoints <- factor(y$samples$timepoint)

logCPMs_corrected <- limma::removeBatchEffect(logCPMs, group = timepoints)
# logCPMs_corrected  <- logCPMs

# Run PCA
pca <- prcomp(t(logCPMs_corrected))

# Combine PCA coordinates with the metadata from the DGEList
to_plot <- data.frame(pca$x, y$samples)

# Calculate how many % of total variance is explained by each principal component
percentVar <- pca$sdev^2 / sum( pca$sdev^2 )*100

# We focus here on PC1 and PC2
use.pcs <- c(1,2)
labs <- paste0(paste0("PC", use.pcs, " - "), paste0("Var.expl = ", round(percentVar[use.pcs], 2), "%"))

ggplot(to_plot, aes(x=PC1, y=PC2,  color=timepoint)) + 
  geom_point(size=3) +
  xlab(labs[1]) + ylab(labs[2])
# ggsave("/mnt/nas-safu/analysis/PhDsdigiove/method_coAcces/data/BALL/communities_broad/AFTER_correction_PCA.png" , height = 5 , width = 9 , units = 'in' , dpi = 300)

logCPMs_corrected_filtered <- logCPMs_corrected #|> as.data.frame()

#logCPMs_corrected_filtered$musmusculus_heart_13_5d_rep2  <- logCPMs_corrected_filtered$musmusculus_heart_13_5d_rep1

cor_matrix <- cor(logCPMs_corrected_filtered , method = "pearson") # Use 'spearman' or 'kendall' if needed


# Plot the heatmap
ph  <- pheatmap(cor_matrix, 
         #display_numbers = TRUE, # Show correlation values on the heatmap
         clustering_method = "complete", # Clustering method
         color = colorRampPalette(c("blue", "white", "red"))(100)) # Color gradient
# png("/mnt/nas-safu/analysis/PhDsdigiove/method_coAcces/data/BALL/communities_broad/AFTER_correction_Correlation.png" , height = 5 , width = 9 , units = 'in' , res = 300)
ph
# dev.off()

logCPMs_corrected_filtered_noextra <- logCPMs_corrected_filtered[!grepl("chrX", rownames(logCPMs_corrected_filtered)) , ]
logCPMs_corrected_filtered_noextra <- logCPMs_corrected_filtered_noextra[!grepl("chrY", rownames(logCPMs_corrected_filtered_noextra)) , ]
logCPMs_corrected_filtered_noextra <- logCPMs_corrected_filtered_noextra[!grepl("gl0", rownames(logCPMs_corrected_filtered_noextra)) , ]
logCPMs_corrected_filtered_noextra <- logCPMs_corrected_filtered_noextra[!grepl("random", rownames(logCPMs_corrected_filtered_noextra)) , ]
logCPMs_corrected_filtered_noextra <- logCPMs_corrected_filtered_noextra[!grepl("chrM", rownames(logCPMs_corrected_filtered_noextra)) , ]

row_mean = rowMeans2(logCPMs_corrected_filtered_noextra)
row_var = rowVars(logCPMs_corrected_filtered_noextra)

data.frame("rowMean" = row_mean , "rowVar" = row_var ) |> 
  ggplot(aes(x = rowMean , y = rowVar))+
  geom_point() +
  geom_hline(yintercept = 1.0) +
  theme_classic()
# ggsave("/mnt/nas-safu/analysis/PhDsdigiove/method_coAcces/data/BALL/communities_broad/Means_vs_Variance_BALL_nostrangeCHR.png" , height = 5 , width = 9 , units = 'in' , dpi = 300)


cor_matrix <- cor(logCPMs_corrected_filtered_noextra , method = "pearson") # Use 'spearman' or 'kendall' if needed

# Plot the heatmap
# png("/mnt/nas-safu/analysis/PhDsdigiove/method_coAcces/pictures/cell_report/batch_correction/AFTER_correction_PCA.png" , height = 5 , width = 9 , units = 'in' , res = 300)
out  <- pheatmap(cor_matrix, 
         #display_numbers = TRUE, # Show correlation values on the heatmap
         clustering_method = "complete", # Clustering method
         color = colorRampPalette(c("blue", "white", "red"))(100)) # Color gradient
# dev.off()

logCPMs_var  <- logCPMs_corrected_filtered_noextra[row_var > 1.0,] |> as.data.frame()
# logCPMs_var |> rownames_to_column('peaks')  |> write_delim( "/mnt/nas-safu/analysis/PhDsdigiove/method_coAcces/data/BALL/allPatientsPatientsNOMERGED_Ball_Multicov_nfcore.tsv" , delim = "\t")
