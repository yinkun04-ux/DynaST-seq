library(Seurat)

filter_rna_data <- function(rna_data, total_cells, expression_threshold = 0.03, cell_threshold = 0.01) {
        total_expression <- colSums(rna_data)
        cells_with_expression <- colSums(rna_data > 0)
        filter_gene <- (total_expression > (total_cells * expression_threshold)) & 
                        (cells_with_expression > (total_cells * cell_threshold))
        return(rna_data[filter_gene, ])
}

seurat_data <- readRDS("/home/disk/chengrui/workspace/spatial_transcriptomics/DynamicST/final_summary_script/result/4_multisample_integration/ureter_added/results/all_RNA.rds")

for (sample_name in unique(seurat_data$orig.ident)) {
        
        region_data <- subset(seurat_data, orig.ident == sample_name)
        
        newRNA_data <- GetAssayData(region_data, assay = "NewRNA", slot = "counts")
        oldRNA_data <- GetAssayData(region_data, assay = "OldRNA", slot = "counts")
        totalRNA_data <- newRNA_data + oldRNA_data
        
        total_cells <- ncol(totalRNA_data)
        filtered_totalRNA_data <- filter_rna_data(totalRNA_data, total_cells)
        filtered_newRNA_data <- filter_rna_data(newRNA_data, total_cells)
        filtered_oldRNA_data <- filter_rna_data(oldRNA_data, total_cells)
        
        write.csv(t(as.matrix(filtered_totalRNA_data)),
                file = paste0("/home/disk/chengrui/workspace/spatial_transcriptomics/DynamicST/final_summary_script/6_transcription_factors_gy/data/raw_counts/",
                                sample_name, "_total.csv"))
        write.csv(t(as.matrix(filtered_newRNA_data)),
                file = paste0("/home/disk/chengrui/workspace/spatial_transcriptomics/DynamicST/final_summary_script/6_transcription_factors_gy/data/raw_counts/",
                                sample_name, "_new.csv"))
        write.csv(t(as.matrix(filtered_oldRNA_data)),
                file = paste0("/home/disk/chengrui/workspace/spatial_transcriptomics/DynamicST/final_summary_script/6_transcription_factors_gy/data/raw_counts/",
                                sample_name, "_old.csv"))
}
