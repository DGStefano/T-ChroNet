library(TCseq)
library(dplyr)
library(cluster)

# --- ARGUMENT PARSING ---
args <- commandArgs(trailingOnly = TRUE)
bed_path <- as.numeric(args[1])

# --- REMOVE GARBAGE ---
gc(full = TRUE)

# --- CONFIGURATION ---
BED_DIR <- bed_path
BAM_DIR <- "/mnt/nas-safu01/analysis/PhDsdigiove/method_coAcces/data/BALL/preliminary_analysis/results_atac/bwa/merged_library/"

data_tcseq <- readr::read_delim(BED_DIR) |> as.data.frame()
# --- DATA PREPARATION ---
gf_full <- peakreference(data = data_tcseq, pattern = "bed")

bamfiles <- data.frame(
    sampleid = c("patient_4", "patient_5", "patient_6", "patient_7", "patient_8", "patient_9", "patient_7", "patient_10", "patient_11", "patient_12", "patient_13", "patient_14", "patient_16", "patient_17", "patient_26", "patient_27", "patient_19", "patient_18", "patient_20", "patient_21", "patient_22", "patient_23", "patient_24", "patient_25", "patient_26", "patient_14", "patient_27", "patient_28", "patient_29", "patient_30", "patient_31", "patient_32", "patient_33"),
    timepoint = c("H","H","H","H","H","H","PRI","PRI","PRI","PRI","PRI","PRI","PRI","PRI","PRI","PRI","PRI","REM","REM","REM","REM","REM","REM","REM","REM","REL","REL","REL","REL","REL","REL","REL","REL"),
    group = c("H","H","H","H","H","H","PRI","PRI","PRI","PRI","PRI","PRI","PRI","PRI","PRI","PRI","PRI","REM","REM","REM","REM","REM","REM","REM","REM","REL","REL","REL","REL","REL","REL","REL","REL"),
    bamfile = c("patient_4_Healthy_ATACseq_REP1.mLb.clN.sorted.bam", "patient_5_Healthy_ATACseq_REP1.mLb.clN.sorted.bam", "patient_6_Healthy_ATACseq_REP1.mLb.clN.sorted.bam", "patient_7_Healthy_ATACseq_REP1.mLb.clN.sorted.bam", "patient_8_Healthy_ATACseq_REP1.mLb.clN.sorted.bam", "patient_9_Healthy_ATACseq_REP1.mLb.clN.sorted.bam", "patient_7_Primary_ATACseq_REP1.mLb.clN.sorted.bam", "patient_10_Primary_ATACseq_REP1.mLb.clN.sorted.bam", "patient_11_Primary_ATACseq_REP1.mLb.clN.sorted.bam", "patient_12_Primary_ATACseq_REP1.mLb.clN.sorted.bam", "patient_13_Primary_ATACseq_REP1.mLb.clN.sorted.bam", "patient_14_Primary_ATACseq_REP1.mLb.clN.sorted.bam", "patient_16_Primary_ATACseq_REP1.mLb.clN.sorted.bam", "patient_17_Primary_ATACseq_REP1.mLb.clN.sorted.bam", "patient_26_Primary_ATACseq_REP1.mLb.clN.sorted.bam", "patient_27_Primary_ATACseq_REP1.mLb.clN.sorted.bam", "patient_19_Primary_ATACseq_REP1.mLb.clN.sorted.bam", "patient_18_Remission_ATACseq_REP1.mLb.clN.sorted.bam", "patient_20_Remission_ATACseq_REP1.mLb.clN.sorted.bam", "patient_21_Remission_ATACseq_REP1.mLb.clN.sorted.bam", "patient_22_Remission_ATACseq_REP1.mLb.clN.sorted.bam", "patient_23_Remission_ATACseq_REP1.mLb.clN.sorted.bam", "patient_24_Remission_ATACseq_REP1.mLb.clN.sorted.bam", "patient_25_Remission_ATACseq_REP1.mLb.clN.sorted.bam", "patient_26_Remission_ATACseq_REP1.mLb.clN.sorted.bam", "patient_14_Relapse_ATACseq_REP1.mLb.clN.sorted.bam", "patient_27_Relapse_ATACseq_REP1.mLb.clN.sorted.bam", "patient_28_Relapse_ATACseq_REP1.mLb.clN.sorted.bam", "patient_29_Relapse_ATACseq_REP1.mLb.clN.sorted.bam", "patient_30_Relapse_ATACseq_REP1.mLb.clN.sorted.bam", "patient_31_Relapse_ATACseq_REP1.mLb.clN.sorted.bam", "patient_32_Relapse_ATACseq_REP1.mLb.clN.sorted.bam", "patient_33_Relapse_ATACseq_REP1.mLb.clN.sorted.bam")
)

# --- EXECUTION STEPS (Timed by the Python Master script) ---

tca <- TCA(design = bamfiles, genomicFeature = gf_full)
tca <- countReads(tca, dir = BAM_DIR)
tca <- DBanalysis(tca, filter.type = 'raw', filter.value = 5)

tca <- timecourseTable(tca , norm.method = 'cpm' , filter = FALSE)

t <- tcTable(tca)
# Find best cluster through silhouette score
res_k <- fviz_nbclust(std_data, kmeans, method = "silhouette", k.max = 10)

# Extract the k with the highest average silhouette width
best_k <- as.numeric(res_k$data$clusters[which.max(res_k$data$y)])

# Select the optimal k (first max or Tibs2001 method)
best_k <- maxSE(gap_stat$Tab[, "gap"], gap_stat$Tab[, "SE.sim"], method = "Tibs2001SEmax")

# Final clustering with the selected best_k
tca <- timeclust(tca, algo = 'km', k = best_k, standardize = TRUE)