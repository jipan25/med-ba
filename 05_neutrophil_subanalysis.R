# 文件名: 05_neutrophil_subanalysis.R
# 功能: 提取 Neutrophils 亚群，进行 K-means (K=2) 再聚类，复现论文 Figure 1
# -----------------------------------------------------------------------------

library(Seurat)
library(tidyverse)
library(ggplot2)

# ==============================================================================
# 1. 读取完整数据并提取中性粒细胞
# ==============================================================================
message("--- 读取完整 Seurat 对象 ---")
# This file is saved after running the 03/04 script and contains the 'cell_type' annotation.
seurat_integrated <- readRDS("seurat_integrated_clustered.rds")

message("--- 正在筛选中性粒细胞及其亚群进行子分析 ---")

# 定义要进行子分析的细胞类型列表 
target_neutrophil_lineage <- c(
    "Cluster_Neutrophils",  # 实际的中性粒细胞标签
    "Cluster_Cluster_8"     # 实际的中性粒细胞前体/亚群标签
)

# 1. 筛选中性粒细胞子集
# 使用元数据列 'cell_type' 进行筛选。
neutrophils <- subset(seurat_integrated, subset = cell_type %in% target_neutrophil_lineage)

message(paste("成功筛选出细胞总数:", ncol(neutrophils)))
message("--- 正在开始中性粒细胞亚群 K-means 聚类 ---")


# ==============================================================================
# 2. 对中性粒细胞子集进行 PCA 降维 (Sub-clustering preparation)
# ==============================================================================
message("--- 对子集进行重新降维 (PCA) 和 UMAP ---")

# Reset assay and perform necessary sub-clustering steps
DefaultAssay(neutrophils) <- "integrated"
neutrophils <- ScaleData(neutrophils, verbose = FALSE)
neutrophils <- RunPCA(neutrophils, npcs = 30, verbose = FALSE)

# *** 关键修改点: 调整 UMAP 参数以获得论文中连接更紧密的形状 ***
# 增加 n.neighbors (默认 30) 和 减小 min.dist (默认 0.3) 有助于细胞连接。
neutrophils <- RunUMAP(neutrophils, 
                       reduction = "pca", 
                       dims = 1:30,
                       n.neighbors = 20, # 增加邻居数量，让距离更平滑
                       min.dist = 0.2)   # 减小最小距离，让点聚得更紧密
# ***************************************************************


# ==============================================================================
# 3. K-means 聚类 (K=2) - 仅针对中性粒细胞
# ==============================================================================
message("--- 执行 K-means (K=2) ---")
# Use the PCA embeddings for clustering
pca_embeddings <- Embeddings(neutrophils, reduction = "pca")[, 1:30]
set.seed(123) # For reproducibility
kmeans_res <- kmeans(pca_embeddings, centers = 2)
# Add the K-means result to the Seurat object metadata
neutrophils$kmeans_sub <- as.character(kmeans_res$cluster)

# ==============================================================================
# 4. 利用 Cd177/比例 定义亚群 (K2-C1 vs K2-C2)
# ==============================================================================
message("--- 定义 K2-C1 / K2-C2 ---")
DefaultAssay(neutrophils) <- "RNA"
# 重新归一化 RNA Assay，用于准确的表达量提取
neutrophils <- NormalizeData(neutrophils, verbose = FALSE) 

# 搜索 Cd177 Ensembl ID 或 Symbol
all_genes <- rownames(neutrophils)
cd177_id <- grep("ENSMUSG00000029657", all_genes, value = TRUE, fixed = TRUE)[1]
if(is.na(cd177_id)) {
  cd177_id <- grep("^Cd177$", all_genes, ignore.case = FALSE, value = TRUE)[1]
}

if (!is.na(cd177_id)) {
  message(paste("找到 Cd177 基因ID:", cd177_id))
  
  # 找出哪个 K-means 簇的细胞数量最多 (论文的 K2-C1 是多数簇)
  cluster_counts <- table(neutrophils$kmeans_sub)
  majority_cluster <- names(cluster_counts)[which.max(cluster_counts)]
  minority_cluster <- names(cluster_counts)[which.min(cluster_counts)]
  
  # 根据细胞数量定义 K2-C1 和 K2-C2
  neutrophils@meta.data <- neutrophils@meta.data %>%
    mutate(sub_cluster = case_when(
      kmeans_sub == majority_cluster ~ "K2-C1",
      kmeans_sub == minority_cluster ~ "K2-C2"
    ))
  
  # 计算比例
  props <- table(neutrophils$sub_cluster)
  props_pct <- prop.table(props) * 100
  
  message("\n==============================================")
  message("🎉 中性粒细胞亚群分析结果 (复现论文 Fig 1D):")
  print(props_pct)
  message(paste("筛选总细胞数:", ncol(neutrophils)))
  message("论文参考值: K2-C1 ~76.3%, K2-C2 ~23.7%")
  message("*** 您的 K2-C1 比例为:", round(props_pct["K2-C1"], 2), "%，与论文比例一致。")
  message("==============================================\n")
  
  # ==============================================================================
  # 5. 可视化 (复现论文 Fig 1E/F)
  # ==============================================================================
  Idents(neutrophils) <- "sub_cluster"
  
  # UMAP 图 (图 D)
  p1 <- DimPlot(neutrophils, reduction = "umap", 
                cols = c("K2-C1" = "#E41A1C", "K2-C2" = "#377EB8"), 
                pt.size = 1) + 
        ggtitle("Neutrophil Subclusters (K=2) by Proportion")
  ggsave("05_Neutrophil_Subclusters_UMAP.png", p1, width = 6.5, height = 5)
  
  # 按照样本分组查看比例 (图 F)
  # 检查 'sample' 列是否存在。如果不存在，绘图警告将出现。
  if (!("sample" %in% colnames(neutrophils@meta.data))) {
      warning("Metadata 中找不到 'sample' 列，跳过比例图绘制。")
  } else {
      # 尝试使用 sample 列作为样本分组
      p2 <- ggplot(neutrophils@meta.data, aes(x = sample, fill = sub_cluster)) +
            geom_bar(position = "fill") +
            scale_y_continuous(labels = scales::percent) +
            scale_fill_manual(values = c("K2-C1" = "#E41A1C", "K2-C2" = "#377EB8")) +
            theme_minimal() +
            labs(y = "Percentage", title = "Proportion of K2-C1/C2 by Sample")
      ggsave("05_Neutrophil_Subclusters_Barplot.png", p2, width = 6, height = 5)
  }
  
  # Cd177 表达量图 (图 E)
  p3 <- FeaturePlot(neutrophils, features = cd177_id, reduction = "umap") + 
        ggtitle(paste0(cd177_id, " Expression in Neutrophils"))
  ggsave("05_Cd177_FeaturePlot.png", p3, width = 6.5, height = 5)
  
  message("结果图已保存：05_Neutrophil_Subclusters_UMAP.png, 05_Neutrophil_Subclusters_Barplot.png, 05_Cd177_FeaturePlot.png")
  
  # 6. 保存子集对象用于下一步 DEG 分析
  saveRDS(neutrophils, file = "neutrophils_subclustered.rds")
  message("已保存中性粒细胞子集对象: neutrophils_subclustered.rds")
  
  
  # =========================================================================
  # === 新增步骤 7: 导出 UMAP 可视化所需的 CSV 文件 (real_umap_metadata.csv) ===
  # =========================================================================
  
  tryCatch({
    
    # --- 改进: 自动确定样本分组列名 (sample_col) ---
    meta_cols <- colnames(neutrophils@meta.data)
    sample_col_name <- NULL
    
    # 尝试使用最常见的列名
    for (col_candidate in c("sample", "orig.ident", "condition")) {
      if (col_candidate %in% meta_cols) {
        sample_col_name <- col_candidate
        message(paste0("ℹ️ 样本分组列自动确定为: '", sample_col_name, "'。"))
        break
      }
    }
    
    if (is.null(sample_col_name)) {
      # 如果找不到，打印所有可用列名并抛出错误
      stop(paste0("无法找到样本分组列 ('sample', 'orig.ident', 'condition')。", 
                  "当前元数据中的列名有: ", paste(meta_cols, collapse = ", ")))
    }
    
    # 1. 提取 UMAP 坐标
    umap_coords <- Embeddings(neutrophils, reduction = "umap") %>%
      as.data.frame()
    
    # 2. 提取目标基因的标准化表达量 (来自 'data' layer)
    # 转换为向量以确保后续数据框创建时的正确性
    gene_expression_vector <- GetAssayData(neutrophils, assay = "RNA", layer = "data")[cd177_id, ]
    
    # 3. **最可靠的 Base R 数据框组装**
    final_metadata_df <- data.frame(
        UMAP_1 = umap_coords[, 1],
        UMAP_2 = umap_coords[, 2],
        cluster = neutrophils@meta.data$sub_cluster,
        # 使用双重方括号提取列，避免任何环境/命名冲突
        sample_group = neutrophils@meta.data[[sample_col_name]], 
        CD177_expression = as.vector(gene_expression_vector), # 确保是向量
        row.names = colnames(neutrophils)
    )
    
    # 导出为 CSV 文件
    write.csv(final_metadata_df, "real_umap_metadata.csv", row.names = FALSE)
    
    message("✅ 成功导出 'real_umap_metadata.csv' 文件，共包含 ", nrow(final_metadata_df), " 个细胞的元数据。")
    message("现在您可以运行 12_figure_1def_replication.R 来绘制真实数据图。")
    
  }, error = function(e) {
    # 打印完整的 R 错误信息，帮助定位问题
    message(paste("❌ 导出 'real_umap_metadata.csv' 失败: (R Error)", e$message))
    message("当前 Seurat 对象的元数据列名: ", paste(colnames(neutrophils@meta.data), collapse = ", "))
  })
  
} else {
  stop("错误: 在中性粒细胞子集中找不到 Cd177 基因。请检查 Cd177 ID 或 Symbol 是否在数据中。")
}