# 文件名: 03_visualization_and_annotation.R (修正版 V2)

library(Seurat)
library(ggplot2)
# 文件名: 03_visualization_and_annotation.R (修正版 V3)

# ... (保持不变)
library(dplyr)
library(cowplot)
library(EnsDb.Mmusculus.v79) # 👈 确保加载
# ... (保持不变)

# ==============================================================================
# 步骤 1: 读取数据和 UMAP 可视化
# ==============================================================================

seurat_integrated <- readRDS("seurat_integrated_clustered.rds")
DefaultAssay(seurat_integrated) <- "RNA" # 切换回 RNA 原始数据

message("--- 正在进行 UMAP 可视化 ---")

# ------------------------------------------------------------------------------
# 关键修正: 合并图层以进行差异表达分析
# ------------------------------------------------------------------------------
message("--- 正在清理并合并 RNA Assay 中的图层 ---")

# 1. 临时保存 RNA 原始 Counts
rna_counts <- GetAssayData(seurat_integrated, assay = "RNA", layer = "counts")

# 2. 临时切换默认 Assay，以便删除 RNA Assay
if (DefaultAssay(seurat_integrated) == "RNA") {
    # 假设 'integrated' Assay 存在
    DefaultAssay(seurat_integrated) <- "integrated" 
}

# 3. 移除旧的 RNA Assay 并重新添加
seurat_integrated[["RNA"]] <- NULL
seurat_integrated[["RNA"]] <- CreateAssayObject(counts = rna_counts)

# 4. 将默认 Assay 切换回 RNA
DefaultAssay(seurat_integrated) <- "RNA" 

# 5. 重新进行标准化，生成 'data' layer
seurat_integrated <- NormalizeData(seurat_integrated, assay = "RNA")

# 6. 执行 JoinLayers：显式指定使用 Seurat 包中的 JoinLayers
#    ----------------------------------------------------------------------
# seurat_integrated <- JoinLayers(object = seurat_integrated, assay = "RNA") 
#    ----------------------------------------------------------------------


# ==============================================================================
# 步骤 2: 寻找 Marker 基因
# ==============================================================================

message("--- 正在寻找所有聚类簇的 Marker 基因 (FindAllMarkers) ---")

# 使用 FindAllMarkers 查找每个聚类簇的显著高表达基因
all_markers <- FindAllMarkers(seurat_integrated, 
                              only.pos = TRUE, 
                              min.pct = 0.25, 
                              logfc.threshold = 0.25)

# ... (后面的保存和绘制部分不变) ...

# 保存所有的 Marker 基因表
write.csv(all_markers, file = "03_All_Cluster_Markers.csv", row.names = FALSE)
message("所有 Marker 基因已保存为 03_All_Cluster_Markers.csv")

# 提取每个聚类簇的前 5 个 top marker 基因
top5_markers <- all_markers %>%
    group_by(cluster) %>%
    slice_max(n = 5, order_by = avg_log2FC)

write.csv(top5_markers, file = "03_Top5_Cluster_Markers.csv", row.names = FALSE)
message("每个聚类簇的前 5 个 Marker 基因已保存为 03_Top5_Cluster_Markers.csv")

# ... (步骤 2: FindAllMarkers 之后)

# ... (步骤 2: FindAllMarkers 之后)

# ==============================================================================
# 步骤 3: 绘制关键 Marker 基因的表达图 (使用 Ensembl ID 绘图)
# ==============================================================================

message("--- 正在进行 Marker 基因 Ensembl ID 转换 ---")

# A. ID 转换函数（保持不变，用于获取 Gene Symbol 信息）
convert_id_to_symbol <- function(ensembl_ids) {
    # 移除版本号（例如 .4）
    clean_ids <- gsub("\\.\\d+$", "", ensembl_ids)
    
    symbols <- AnnotationDbi::select(EnsDb.Mmusculus.v79, 
                      keys = clean_ids, 
                      keytype = "GENEID", 
                      columns = "SYMBOL")
    
    converted_list <- setNames(symbols$SYMBOL, symbols$GENEID)
    return(converted_list)
}

# 1. 提取 Cluster 0 的 Top 5 Ensembl ID (用于绘图)
top5_ens_id <- top5_markers %>% 
    dplyr::filter(cluster == 0) %>%
    dplyr::pull(gene)

# ==============================================================================
# 优化后：根据 Gene Symbol 动态查找数据中的 Ensembl ID
# ==============================================================================

# 1. 定义你想要可视化的目标基因 (人类可读的 Symbol)
#    这是你唯一需要修改的“配置列表”
target_symbols <- c("Epcam", "Krt7", "Alb", "Ly6g", "S100a8") 

message("--- 正在将目标 Symbol 映射回数据集中的 Ensembl ID ---")

# 2. 构建全基因组映射字典 (Symbol -> Ensembl ID)
#    (利用已加载的 EnsDb 数据库)
all_ens_ids <- rownames(seurat_integrated) # 获取数据中实际存在的 ID (带版本号)
clean_ens_ids <- gsub("\\.\\d+$", "", all_ens_ids) # 去掉版本号用于查询

# 从数据库批量查询所有 ID 对应的 Symbol
gene_map <- AnnotationDbi::select(EnsDb.Mmusculus.v79, 
                                  keys = clean_ens_ids, 
                                  keytype = "GENEID", 
                                  columns = "SYMBOL")

# 创建一个查找向量: 名字是 Symbol, 值是原始 Ensembl ID (带版本号)
# 我们需要把数据库查到的 Symbol 和数据行名对应起来
# 注意：match 的顺序很重要
match_idx <- match(clean_ens_ids, gene_map$GENEID)
current_symbols <- gene_map$SYMBOL[match_idx]

# 处理可能的 NA (数据库里没查到的基因)
valid_idx <- !is.na(current_symbols)
symbol_to_id_map <- setNames(all_ens_ids[valid_idx], current_symbols[valid_idx])

# 3. 查找目标基因的 ID
#    intersect 确保只查找字典里有的基因
found_symbols <- intersect(target_symbols, names(symbol_to_id_map))
key_ens_id_literature <- symbol_to_id_map[found_symbols]

# 4. 报告查找结果
missing_symbols <- setdiff(target_symbols, found_symbols)
if(length(missing_symbols) > 0) {
  warning(paste("以下基因未在数据集中找到 (可能被过滤或名称不同):", paste(missing_symbols, collapse=", ")))
}

message("成功找到以下 Marker ID:")
print(key_ens_id_literature)

# ... (接下来的 key_markers_plotting 合并逻辑保持不变)

# 3. 最终绘图列表：使用 Ensembl ID 确保 FeaturePlot 找到特征
key_markers_plotting <- unique(c(top5_ens_id, key_ens_id_literature))
# 确保只保留数据中存在的 ID
key_markers_plotting <- intersect(key_markers_plotting, rownames(seurat_integrated))

# 4. 转换 ID 为 Symbol，用于打印和注释 (关键步骤)
message("\n--- Cluster 0 Top 5 Marker (Gene Symbol) ---")
# 仅将 Top 5 Ensembl ID 转换为可读的 Symbol 并打印
top5_gene_symbol <- convert_id_to_symbol(top5_ens_id)
# 打印结果：
print(top5_gene_symbol)


message("--- 正在绘制关键 Marker 基因表达图 (使用 Ensembl ID) ---")

# 绘制 FeaturePlot
marker_plots <- FeaturePlot(seurat_integrated, features = key_markers_plotting, ncol = 3)
ggsave("03_Key_Marker_FeaturePlots.png", plot = marker_plots, width = 15, height = 10)
message("关键 Marker 基因表达图已保存为 03_Key_Marker_FeaturePlots.png")

message("--- 脚本运行完成 ---")

# ==============================================================================
# 步骤 4: 细胞类型注释 (Cell Type Annotation)
# ==============================================================================
message("--- 正在进行细胞类型注释 ---")

# 1. 定义注释字典 (配置部分)
cluster_annotation_map <- c(
    "0" = "Neutrophils",
    "1" = "Monocytes/Macrophages",
    "2" = "Epithelial/iBEC_Progenitor",
    "3" = "T_Cells",
    "4" = "B_Cells"
)

# 2. **关键修正：使用 Idents() 获取当前聚类，并使用名字进行映射**

# 获取当前细胞的聚类身份 (Idents)
current_idents <- Idents(seurat_integrated) 

# 转换为字符类型以便映射
cluster_names <- as.character(current_idents)

# 创建新标签向量，并保持细胞条形码作为名字
new_labels_named <- setNames(cluster_annotation_map[cluster_names], names(current_idents))


# 3. 处理未定义的 Cluster (健壮性处理)

unmapped_idx <- is.na(new_labels_named)

if (any(unmapped_idx)) {
    # 对于没定义的，生成 "Cluster_X" 格式 (这里的名字是 Cluster ID)
    unmapped_clusters <- unique(cluster_names[unmapped_idx])
    
    # 打印提示信息
    message(paste("提示: 以下 Cluster 未定义具体名称，将保留默认编号:", 
                  paste(unmapped_clusters, collapse = ", ")))
                  
    # 将未映射的标签设为 "Cluster_X"
    new_labels_named[unmapped_idx] <- paste0("Cluster_", cluster_names[unmapped_idx])
}


# 4. **核心修正：使用 AddMetaData 赋值回 Seurat 对象**

# 创建一个 metadata 列
metadata_to_add <- data.frame(cell_type = new_labels_named)
rownames(metadata_to_add) <- names(new_labels_named) # 确保行名是细胞条形码

# 使用 AddMetaData 函数安全地添加新列
seurat_integrated <- AddMetaData(object = seurat_integrated, metadata = metadata_to_add, col.name = "cell_type")


# 5. 更新 Seurat 的默认身份 (Idents) 为 cell_type

# 转换为因子 (Factor) 并指定顺序
defined_levels <- unique(cluster_annotation_map)
all_levels <- unique(c(defined_levels, unique(new_labels_named)))
seurat_integrated$cell_type <- factor(seurat_integrated$cell_type, levels = all_levels)

# 设置新的 Idents
Idents(seurat_integrated) <- "cell_type"

# 保存
saveRDS(seurat_integrated, file = "seurat_integrated_clustered.rds")
message("带注释的 Seurat 对象已保存。")


rm(list = ls())