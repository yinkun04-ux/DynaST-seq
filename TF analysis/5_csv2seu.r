library(dplyr)
library(Seurat)
library(tidyverse)
library(future)
options(future.global.maxSize = 1200000000)
library(pheatmap)

process_auc_matrix <- function(path, sample_name) {
  AUCmatrix <- read.csv(paste0(path, '/', sample_name, '_SCENIC.csv'), row.names = 1, check.names = FALSE)
  AUCmatrix <- t(AUCmatrix)
  rownames(AUCmatrix) <- gsub("\\(\\+\\)$", "", rownames(AUCmatrix))
  return(scale(AUCmatrix))
}

prepare_tf_data <- function(matrix) {
  x <- strsplit(row.names(matrix), "_")
  TFs <- sapply(x, function(y) y[1])
  matrix$TF <- TFs
  matrix <- matrix[!duplicated(matrix[, ncol(matrix)]), ]
  row.names(matrix) <- matrix$TF
  matrix1 <- matrix[, -ncol(matrix)]
  return(matrix1)
}

calculate_mean_expression <- function(split_list, matrix1, timepoints = c("0h", "4h", "12h", "24h")) {
  means <- sapply(timepoints, function(time) {
    expr_matrix <- as.data.frame(split_list[[time]]@assays[["RNA"]]$data)[row.names(matrix1),]
    apply(expr_matrix, 1, mean)
  })
  return(data.frame(means, row.names = row.names(matrix1)))
}

create_seurat_object <- function(matrix) {
  seurat_obj <- CreateSeuratObject(matrix, project = "scenic")
  seurat_obj[['RNA']]$data <- seurat_obj[['RNA']]$counts
  return(seurat_obj)
}

# Load AUC matrix and Seurat object
args <- commandArgs(T)
sample_name1 <- args[1]
path <- "/home/disk/chengrui/workspace/spatial_transcriptomics/DynamicST/final_summary_script/6_transcription_factors_gy/result/pic"
AUCmatrix <- process_auc_matrix(path, sample_name1)
seurat_obj <- create_seurat_object(AUCmatrix)

meta_data <- read.csv("/home/disk/chengrui/workspace/spatial_transcriptomics/DynamicST/final_summary_script/result/4_multisample_integration/ureter_added/meta_data.csv", row.names = 1)
meta_data <- meta_data[rownames(seurat_obj@meta.data),]
seurat_obj@meta.data <- cbind(seurat_obj@meta.data, meta_data[, "time"])
colnames(seurat_obj@meta.data) <- c(colnames(seurat_obj@meta.data)[1:(ncol(seurat_obj@meta.data) - 1)], "time")
seurat_obj$time <- paste0(seurat_obj$time, "h")
seurat_obj$time <- as.character(seurat_obj$time)

split.list <- SplitObject(seurat_obj, split.by = "time")
saveRDS(split.list, file = paste0(path, '/', sample_name1, '_split.rds'))

get_markers <- function(seurat_obj, timepoints, threshold = 0.58) {
    markers_list <- list()
    for (timepoint in timepoints) {
        markers_list[[timepoint]] <- FindMarkers(seurat_obj, ident.1 = timepoint, ident.2 = '0h', min.pct = 0.1, logfc.threshold = threshold, test.use = "wilcox")
    }
    return(markers_list)
}

timepoints <- c('4h', '12h', '24h')
markers <- get_markers(seurat_obj, timepoints)

# Filter markers by log fold change threshold
filter_markers <- function(marker_list, logfc_threshold = 0.58) {
    pos_markers <- lapply(marker_list, function(x) subset(x, avg_log2FC > logfc_threshold))
    neg_markers <- lapply(marker_list, function(x) subset(x, avg_log2FC < -logfc_threshold))
    
    pos_table <- do.call(union, lapply(pos_markers, rownames))
    neg_table <- do.call(union, lapply(neg_markers, rownames))
    
    pos_table <- as.data.frame(pos_table)
    neg_table <- as.data.frame(neg_table)
    colnames(pos_table) <- "marker"
    colnames(neg_table) <- "marker"
    
    return(list(pos = pos_table, neg = neg_table))
}

filtered_markers <- filter_markers(markers)

# Calculate mean expression for positive and negative markers
all_merge_expr_mean_pos <- calculate_mean_expression(split.list, AUCmatrix, timepoints)
all_merge_expr_mean_neg <- calculate_mean_expression(split.list, AUCmatrix, timepoints)

# Process for heatmap visualization
prepare_heatmap_data <- function(all_merge_expr_mean) {
    all_merge_expr_mean$TFname <- row.names(all_merge_expr_mean)
    mat <- all_merge_expr_mean[!duplicated(all_merge_expr_mean[, "TFname"]), ]
    row.names(mat) <- mat$TFname
    mat1 <- mat[, -ncol(mat)]
    return(mat1)
}

# Prepare positive and negative marker data
heatmap_data_pos <- prepare_heatmap_data(all_merge_expr_mean_pos)
heatmap_data_neg <- prepare_heatmap_data(all_merge_expr_mean_neg)

# Save the data for further analysis
write.table(heatmap_data_pos, file = paste0(path, '/', sample_name1, '_DE1.5posforpheap.txt'), sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
write.table(heatmap_data_neg, file = paste0(path, '/', sample_name1, '_DE1.5negforpheap.txt'), sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)

# Plot heatmap
pheatmap(heatmap_data_pos, scale = "row", border = FALSE, fontsize_col = 11, cluster_rows = FALSE, cluster_cols = FALSE, fontsize = 11, cellwidth = 30, cellheight = 8, filename = paste0(path, '/', sample_name1, '_targetgenemeanfortime_pos.pdf'))
pheatmap(heatmap_data_neg, scale = "row", border = FALSE, fontsize_col = 11, cluster_rows = FALSE, cluster_cols = FALSE, fontsize = 11, cellwidth = 30, cellheight = 8, filename = paste0(path, '/', sample_name1, '_targetgenemeanfortime_neg.pdf'))

# Save mean expression data for time points
write.table(all_merge_expr_mean_pos, file = paste0(path, '/', sample_name1, '_meanfortime_pos.txt'), sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
write.table(all_merge_expr_mean_neg, file = paste0(path, '/', sample_name1, '_meanfortime_neg.txt'), sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
