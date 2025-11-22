# 文件名: 04_differential_expression_and_enrichment.R

# ==============================================================================
# 1. 环境准备和数据加载
# ==============================================================================
message("--- 正在加载环境和数据 ---")
library(Seurat)
library(dplyr)
library(ggplot2)

# 确保安装并加载富集分析包 (只需要运行一次)
# if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
# BiocManager::install(c("clusterProfiler", "org.Mm.eg.db"), dependencies = TRUE)

library(clusterProfiler)
library(org.Mm.eg.db) # 小鼠基因组注释包

# 加载已聚类和注释的 Seurat 对象
seurat_integrated <- readRDS("seurat_integrated_clustered.rds")


# ==============================================================================
# 2. 差异表达分析 (Differential Expression, DE)
# ==============================================================================

# 目标：比较 Cluster 0 (Neutrophils) 中，两个样本之间的基因表达差异。

# 检查样本条件列名 (假设是 orig.ident)
message("检查 Seurat 对象中的样本/条件信息...")
print(table(seurat_integrated$orig.ident, seurat_integrated$cell_type))

# 设置比较的身份 (Idents) 为样本 (orig.ident)
Idents(seurat_integrated) <- "orig.ident"

# 提取 Cluster 0 的细胞 (根据您的日志，Neutrophils 细胞类型位于 cell_type 列)
neutrophil_cells <- subset(seurat_integrated, subset = cell_type == "Neutrophils")

# ⚠️ 关键步骤：使用 FindMarkers 比较两个样本的差异
# 修正：将 'sample2'/'sample1' 替换为实际的样本标识符 'RRV_5d' 和 'NC_5d'

message("\n--- 正在对 Cluster 0 (Neutrophils) 进行差异表达分析 (RRV_5d vs NC_5d) ---")
de_markers_c0 <- FindMarkers(neutrophil_cells, 
                             ident.1 = "RRV_5d", # 处理组
                             ident.2 = "NC_5d",  # 对照组
                             min.pct = 0.05,        # 降低最小表达比例要求 (如 5%)
                             logfc.threshold = 0.0, # 确保 FindMarkers 包含所有有差异的基因
                             verbose = FALSE)

# 清理结果并保存
de_markers_c0 <- de_markers_c0 %>% 
  tibble::rownames_to_column("gene") %>%
  dplyr::arrange(desc(avg_log2FC))

write.csv(de_markers_c0, "04_Cluster0_DE_Markers_RRV_5d_vs_NC_5d.csv", row.names = FALSE)
message("Cluster 0 的差异表达基因已保存为 04_Cluster0_DE_Markers_RRV_5d_vs_NC_5d.csv")


# ==============================================================================
# 3. 富集分析 (Enrichment Analysis)
# ==============================================================================

# 检查 FindMarkers 结果是否为空
if (nrow(de_markers_c0) == 0) {
  message("\n⚠️ 警告: FindMarkers 返回空结果，无法进行差异表达和富集分析。请尝试放宽 min.pct 或 logfc.threshold 阈值。")
} else {
  
  # A. 提取上调基因 (LogFC > 0.5 且 p_val_adj < 0.05)
  upregulated_genes <- de_markers_c0 %>%
    dplyr::filter(avg_log2FC > 0.5 & p_val_adj < 0.05) %>%
    head(500) # 取 Top 500 个最显著的上调基因

  if(nrow(upregulated_genes) == 0) {
    message("\n⚠️ 警告: 没有足够的显著上调基因进行富集分析 (LogFC > 0.5 & p_val_adj < 0.05)。请调整筛选阈值。")
  } else {
    # B. 将 Ensembl ID 转换为 ENTREZ ID (富集分析需要)
    ensembl_ids <- gsub("\\.\\d+$", "", upregulated_genes$gene)
    
    entrez_ids <- AnnotationDbi::mapIds(org.Mm.eg.db, 
                                        keys = ensembl_ids, 
                                        column = "ENTREZID", 
                                        keytype = "ENSEMBL", 
                                        multiVals = "first")
    
    entrez_ids <- na.omit(entrez_ids)
    
    message(paste0("\n--- 正在对 ", length(entrez_ids), " 个上调基因进行 GO/KEGG 富集分析 ---"))
    
    # C. GO (Gene Ontology) 富集分析
    go_results <- enrichGO(gene = entrez_ids,
                          OrgDb = org.Mm.eg.db,
                          keyType = "ENTREZID",
                          ont = "BP", # Biological Process
                          pAdjustMethod = "BH",
                          qvalueCutoff = 0.05)
    
    # D. KEGG 富集分析
    kegg_results <- enrichKEGG(gene = entrez_ids,
                              organism = 'mmu', # 小鼠
                              pAdjustMethod = "BH",
                              qvalueCutoff = 0.05)
    
    # E. 保存和可视化结果
    write.csv(as.data.frame(go_results), "04_Cluster0_GO_Enrichment.csv", row.names = FALSE)
    write.csv(as.data.frame(kegg_results), "04_Cluster0_KEGG_Enrichment.csv", row.names = FALSE)
    message("GO/KEGG 富集结果已保存。")
    
    # 可视化 Top 10 通路
    if (!is.null(go_results)) {
      go_plot <- barplot(go_results, showCategory = 10, title = "GO Biological Process Enrichment")
      ggsave("04_Cluster0_GO_Enrichment_Plot.png", plot = go_plot, width = 10, height = 7)
      message("GO 富集图已保存。")
    }
  }
}

message("\n--- 04_differential_expression_and_enrichment.R 脚本运行完成 ---")