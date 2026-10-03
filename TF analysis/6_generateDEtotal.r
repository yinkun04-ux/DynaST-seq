library(Matrix)
library(dplyr)
library(patchwork)
library(ggplot2)
library(RColorBrewer)
library(cowplot)
library(pheatmap)

prepare_data <- function(sample_name, file_suffix) {
    read.delim(paste0(res_path, sample_name, file_suffix))
}

process_matrix <- function(new_matrix, old_matrix, common_rows) {
    new_matrix <- new_matrix[common_rows,]
    old_matrix <- old_matrix[common_rows,]
    return(list(new_matrix = new_matrix, old_matrix = old_matrix))
}

create_heatmap <- function(matrix_data, annotation_data, filename) {
    pheatmap(matrix_data, scale = "row", annotation_col = annotation_data, border = FALSE,
            fontsize_col = 11, cluster_cols = FALSE, fontsize = 11, cellwidth = 30, cellheight = 8, 
            filename = filename)
}

z_score_normalize <- function(x) {
    return((x - mean(x)) / sd(x))
}

args <- commandArgs(T)
sample_name <- args[1]
res_path = "/home/disk/chengrui/workspace/时空转录组/DynamicST/最终需求的汇总脚本/6_转录因子gy/result/pic/"

pheapnewmatrixpos <- prepare_data(sample_name, '_new_DE1.5posforpheap.txt')
pheapnewmatrixneg <- prepare_data(sample_name, '_new_DE1.5negforpheap.txt')
pheapoldmatrixpos <- prepare_data(sample_name, '_old_DE1.5posforpheap.txt')
pheapoldmatrixneg <- prepare_data(sample_name, '_old_DE1.5negforpheap.txt')
newmeanmatrix <- prepare_data(sample_name, '_new_meanfortime.txt')
oldmeanmatrix <- prepare_data(sample_name, '_old_meanfortime.txt')

allpos <- union(rownames(pheapnewmatrixpos), rownames(pheapoldmatrixpos)) %>% as.data.frame()
allneg <- union(rownames(pheapnewmatrixneg), rownames(pheapoldmatrixneg)) %>% as.data.frame()

# Filter markers for both new and old mean matrices
dat1 <- allpos[allpos[, 1] %in% rownames(newmeanmatrix), ]
dat2 <- dat1[dat1[, 1] %in% rownames(oldmeanmatrix), ]

common_rows <- dat2$dat2
newmeanmatrix1 <- newmeanmatrix[common_rows,]
oldmeanmatrix1 <- oldmeanmatrix[common_rows,]

allmeanmatpos <- cbind(newmeanmatrix1, oldmeanmatrix1)
col <- colnames(allmeanmatpos)
col <- as.data.frame(col)
col$col <- as.factor(col$col)
col$group <- c(rep("new", 4), rep("old", 4))
col$group <- as.factor(col$group)
col <- col[,-1]
row.names(col) <- colnames(allmeanmatpos)
names(col) <- "group"

create_heatmap(allmeanmatpos, col, paste0(res_path, sample_name, "_newoldpos.pdf"))
write.table(allmeanmatpos, file = paste0(res_path, sample_name, "_allmeanmatpos.txt"), sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)

# Process negative markers
dat1 <- allneg[allneg[, 1] %in% rownames(newmeanmatrix), ]
dat2 <- dat1[dat1[, 1] %in% rownames(oldmeanmatrix), ]
newmeanmatrix2 <- newmeanmatrix[dat2$dat2,]
oldmeanmatrix2 <- oldmeanmatrix[dat2$dat2,]
allmeanmatneg <- cbind(newmeanmatrix2, oldmeanmatrix2)

# Create heatmap for negative markers
create_heatmap(allmeanmatneg, col, paste0(res_path, sample_name, "_newoldneg.pdf"))
write.table(allmeanmatneg, file = paste0(res_path, sample_name, "_allmeanmatneg.txt"), sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)

totalmeanmatrix <- newmeanmatrix
allmeanselect <- allmeanmatpos
alltotalnew <- totalmeanmatrix[row.names(allmeanselect),]
alltotalnew <- na.omit(alltotalnew)
write.table(alltotalnew, file = paste0(res_path, sample_name, '_new_allmean.txt'), sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
create_heatmap(alltotalnew, col, paste0(res_path, sample_name, "_new.pdf"))

totalmeanmatrix <- oldmeanmatrix
alltotalold <- totalmeanmatrix[row.names(allmeanselect),]
alltotalold <- na.omit(alltotalold)
write.table(alltotalold, file = paste0(res_path, sample_name, '_old_allmean.txt'), sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
create_heatmap(alltotalold, col, paste0(res_path, sample_name, "_old.pdf"))

new <- alltotalnew
old <- oldmeanmatrix
newtf <- rownames(new)
oldtf <- rownames(old)

# Z-score normalization for comparison
x <- intersect(newtf, oldtf)
if (length(x) > 0) {
    newmat <- new[x,]
    oldmat <- old[x,]

    colnames(newmat) <- paste0("new_", colnames(newmat))
    colnames(oldmat) <- paste0("old_", colnames(oldmat))
    bindmat <- cbind(newmat, oldmat)

    # Z-score normalize
    norm_data <- apply(bindmat, 1, z_score_normalize)
    norm_data <- as.data.frame(t(norm_data))

    annotation_col <- data.frame(Category = factor(rep(c("new", "old"), each = 4)))
    rownames(annotation_col) <- colnames(norm_data)

    annotation_colors <- list(Category = c(new = "#c472a2", old = "#549ad1"))
    
    create_heatmap(norm_data, annotation_col, paste0(res_path, sample_name, "_newold.pdf"))
    write.table(norm_data, paste0(res_path, sample_name, "_newold.txt"), sep = "\t", quote = F)
}
