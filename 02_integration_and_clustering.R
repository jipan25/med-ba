# 文件名: 02_integration_and_clustering.R
# 功能: 对QC后的数据进行标准化、整合 (CCA)、降维 (PCA/UMAP) 和聚类 (Resolution 0.6)。

library(Seurat)
library(tidyverse)

# ==============================================================================
# 步骤 1: 读取数据和标准化
# ==============================================================================

# 读取上一步保存的过滤后的 Seurat 对象
seurat_combined <- readRDS("seurat_combined_filtered.rds")

message(paste("开始处理细胞数:", ncol(seurat_combined)))

# 标准化数据：使用 LogNormalize 方法 (默认)
seurat_list <- SplitObject(seurat_combined, split.by = "sample")
seurat_list <- lapply(seurat_list, function(x) {
    x <- NormalizeData(x)
    # 查找高变基因 (默认 2000 个)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
    return(x)
})

# ... (前面的 NormalizeData 和 FindVariableFeatures 部分不变)

# ==============================================================================
# 步骤 2: 整合 (Canonical Correlation Analysis, CCA)
# ==============================================================================

message("--- 正在使用 CCA 寻找整合锚点 (FindIntegrationAnchors) ---")

# 提取整合特征：合并所有样本的高变基因列表，并确保它是字符向量
integration.features <- SelectIntegrationFeatures(object.list = seurat_list, nfeatures = 2000)
integration.features <- as.character(integration.features) # 关键修正！

# 使用 SelectIntegrationFeatures 提取的基因列表进行整合
anchors <- FindIntegrationAnchors(object.list = seurat_list, 
                                  dims = 1:30, 
                                  anchor.features = integration.features) # 使用修正后的列表

message("--- 正在整合数据 (IntegrateData) ---")

# 将数据整合到一起
seurat_integrated <- IntegrateData(anchorset = anchors, dims = 1:30)

# ... (后面的 RunPCA, RunUMAP, FindClusters 部分不变)

# ==============================================================================
# 步骤 3: 降维和聚类
# ==============================================================================

message("--- 正在进行降维和聚类 ---")

# 切换到整合后的数据 Assays，并进行标准化 (ScaleData)
DefaultAssay(seurat_integrated) <- "integrated"
seurat_integrated <- ScaleData(seurat_integrated, verbose = FALSE)

# 运行 PCA (主成分分析)
seurat_integrated <- RunPCA(seurat_integrated, npcs = 30, verbose = FALSE)

# 运行 UMAP (非线性降维，用于可视化)
seurat_integrated <- RunUMAP(seurat_integrated, reduction = "pca", dims = 1:30)

# 聚类 (Clustering)
# 论文指定：分辨率 0.6
seurat_integrated <- FindNeighbors(seurat_integrated, reduction = "pca", dims = 1:30)
seurat_integrated <- FindClusters(seurat_integrated, resolution = 0.6)

message("--- 聚类完成 ---")
message(paste("最终找到的聚类簇 (Clusters) 数量:", length(unique(seurat_integrated$seurat_clusters))))

# ==============================================================================
# 步骤 4: 保存结果
# ==============================================================================
saveRDS(seurat_integrated, file = "seurat_integrated_clustered.rds")
message("整合后的 Seurat 对象已保存为 seurat_integrated_clustered.rds")

# 打印聚类结果和样本分布 (可选，用于快速检查)
print(table(seurat_integrated$seurat_clusters))
print(table(seurat_integrated$seurat_clusters, seurat_integrated$sample))

rm(list = ls())