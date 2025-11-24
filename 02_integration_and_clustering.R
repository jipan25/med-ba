# 文件名: 02_integration_and_clustering.R (标准版 - 用于全局注释)

library(Seurat)
library(tidyverse)

# 1. 读取数据
seurat_combined <- readRDS("seurat_combined_filtered.rds")
message(paste("处理细胞数:", ncol(seurat_combined)))

# 2. 标准化与整合准备
seurat_list <- SplitObject(seurat_combined, split.by = "sample")
seurat_list <- lapply(seurat_list, function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
    return(x)
})

# 3. CCA 整合
integration.features <- SelectIntegrationFeatures(object.list = seurat_list, nfeatures = 2000)
anchors <- FindIntegrationAnchors(object.list = seurat_list, anchor.features = integration.features)
seurat_integrated <- IntegrateData(anchorset = anchors)

# 4. 降维与聚类 (标准流程)
DefaultAssay(seurat_integrated) <- "integrated"
seurat_integrated <- ScaleData(seurat_integrated, verbose = FALSE)
seurat_integrated <- RunPCA(seurat_integrated, npcs = 30, verbose = FALSE)
seurat_integrated <- RunUMAP(seurat_integrated, reduction = "pca", dims = 1:30)

# 使用标准的分辨率聚类，保留细胞类型的多样性
seurat_integrated <- FindNeighbors(seurat_integrated, reduction = "pca", dims = 1:30)
seurat_integrated <- FindClusters(seurat_integrated, resolution = 0.6)

# 5. 保存
saveRDS(seurat_integrated, file = "seurat_integrated_clustered.rds")
message("标准聚类完成。已保存为 seurat_integrated_clustered.rds")