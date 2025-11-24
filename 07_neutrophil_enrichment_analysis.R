# 文件名: 07_neutrophil_enrichment_analysis.R
# 功能: 对中性粒细胞 K2-C1 vs K2-C2 差异表达基因 (DEG) 进行 GO 和 KEGG 富集分析

library(Seurat)
library(tidyverse)
library(ggplot2)

# ==============================================================================
# 0. 检查和安装富集分析所需包
# ==============================================================================
# 富集分析需要 Bioconductor 包。如果未安装，请在 R 命令行中运行以下代码：
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("clusterProfiler")
# BiocManager::install("org.Mm.eg.db") # 小鼠基因注释包

library(clusterProfiler)
library(org.Mm.eg.db) # 适用于小鼠 (Mus musculus)

# ==============================================================================
# 1. 读取 DEG 结果并筛选
# ==============================================================================
message("--- 读取 DEG 结果 ---")
deg_results <- read.csv("06_Neutrophil_K2_DEG_results.csv")

# 筛选用于富集分析的基因列表
# 使用更严格的阈值来保证富集分析的质量: |logFC| > 0.5 且 p_adj < 0.05
# K2-C1 上调基因 (logFC > 0.5)
deg_up_k2c1 <- deg_results %>%
    filter(avg_log2FC > 0.5, p_val_adj < 0.05)

# K2-C2 上调基因 (logFC < -0.5)
deg_up_k2c2 <- deg_results %>%
    filter(avg_log2FC < -0.5, p_val_adj < 0.05)

message(paste("K2-C1 上调基因数量 (logFC > 0.5):", nrow(deg_up_k2c1)))
message(paste("K2-C2 上调基因数量 (logFC < -0.5):", nrow(deg_up_k2c2)))

# ==============================================================================
# 2. 基因 ID 转换 (Ensembl ID -> Entrez ID)
# clusterProfiler 主要使用 Entrez ID
# ==============================================================================
message("--- 基因 ID 转换: Ensembl ID -> Entrez ID (去除版本号) ---")

# --- K2-C1 上调基因 ---
genes_up_k2c1_raw <- deg_up_k2c1$gene
# FIX: 去除 Ensembl ID 的版本号 (如 .3, .2)
genes_up_k2c1 <- sub("\\..*$", "", genes_up_k2c1_raw)
message(paste("K2-C1 示例基因 (转换前):", paste(head(genes_up_k2c1_raw, 3), collapse = ", ")))
message(paste("K2-C1 示例基因 (转换后):", paste(head(genes_up_k2c1, 3), collapse = ", ")))

eg_up_k2c1 <- bitr(
    geneID = genes_up_k2c1, 
    fromType = "ENSEMBL", 
    toType = c("ENTREZID", "SYMBOL"), 
    OrgDb = org.Mm.eg.db
)

# --- K2-C2 上调基因 ---
genes_up_k2c2_raw <- deg_up_k2c2$gene
# FIX: 去除 Ensembl ID 的版本号 (如 .3, .2)
genes_up_k2c2 <- sub("\\..*$", "", genes_up_k2c2_raw)
message(paste("K2-C2 示例基因 (转换前):", paste(head(genes_up_k2c2_raw, 3), collapse = ", ")))
message(paste("K2-C2 示例基因 (转换后):", paste(head(genes_up_k2c2, 3), collapse = ", ")))

eg_up_k2c2 <- bitr(
    geneID = genes_up_k2c2, 
    fromType = "ENSEMBL", 
    toType = c("ENTREZID", "SYMBOL"), 
    OrgDb = org.Mm.eg.db
)

# 确保至少有 Entrez ID 成功转换
if (nrow(eg_up_k2c1) == 0 && nrow(eg_up_k2c2) == 0) {
    stop("错误: 未能成功将任何基因 ID 转换为 Entrez ID。请检查输入基因ID格式是否正确。")
}

# ==============================================================================
# 3. GO (Gene Ontology) 富集分析
# ==============================================================================
message("--- 执行 GO (BP, MF, CC) 富集分析 ---")

# GO 分析函数
run_go_enrichment <- function(entrez_ids, group_name) {
    message(paste("  -> 正在对", group_name, "进行 GO 分析..."))
    ego <- enrichGO(
        gene = entrez_ids,
        OrgDb = org.Mm.eg.db,
        keyType = "ENTREZID",
        ont = "ALL", # Biological Process, Molecular Function, Cellular Component
        pAdjustMethod = "BH",
        qvalueCutoff = 0.05,
        readable = TRUE # 转换为基因 Symbol
    )
    # 转换为数据框并添加分组信息
    df <- as.data.frame(ego) %>% mutate(Group = group_name)
    # NOTE: clusterProfiler 默认将 GO 本体存储在 'ONTOLOGY' 列
    return(df)
}

# K2-C1 上调基因的 GO 分析
go_k2c1 <- run_go_enrichment(eg_up_k2c1$ENTREZID, "K2-C1_Up")
# K2-C2 上调基因的 GO 分析
go_k2c2 <- run_go_enrichment(eg_up_k2c2$ENTREZID, "K2-C2_Up")

# 合并所有 GO 结果
all_go_results <- bind_rows(go_k2c1, go_k2c2)
write.csv(all_go_results, file = "07_Neutrophil_K2_GO_results.csv", row.names = FALSE)
message("GO 富集结果已保存至: 07_Neutrophil_K2_GO_results.csv")

# ==============================================================================
# 4. KEGG Pathway 富集分析
# ==============================================================================
message("--- 执行 KEGG Pathway 富集分析 ---")

# KEGG 分析函数
run_kegg_enrichment <- function(entrez_ids, group_name) {
    message(paste("  -> 正在对", group_name, "进行 KEGG 分析..."))
    ekegg <- enrichKEGG(
        gene = entrez_ids,
        organism = 'mmu', # 小鼠
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.05
    )
    df <- as.data.frame(ekegg) %>% mutate(Group = group_name)
    return(df)
}

# K2-C1 上调基因的 KEGG 分析
kegg_k2c1 <- run_kegg_enrichment(eg_up_k2c1$ENTREZID, "K2-C1_Up")
# K2-C2 上调基因的 KEGG 分析
kegg_k2c2 <- run_kegg_enrichment(eg_up_k2c2$ENTREZID, "K2-C2_Up")

# 合并所有 KEGG 结果
all_kegg_results <- bind_rows(kegg_k2c1, kegg_k2c2)
write.csv(all_kegg_results, file = "07_Neutrophil_K2_KEGG_results.csv", row.names = FALSE)
message("KEGG 富集结果已保存至: 07_Neutrophil_K2_KEGG_results.csv")

# ==============================================================================
# 5. 可视化 (GO Top 10)
# ==============================================================================
message("--- 生成 GO Top 10 可视化 ---")

# 筛选每个类别和组别的 Top 10
plot_go_data <- all_go_results %>%
    # FIX: 将 ONT 更改为正确的列名 ONTOLOGY
    filter(ONTOLOGY %in% c("BP", "MF", "CC")) %>% 
    # FIX: 将 ONT 更改为正确的列名 ONTOLOGY
    group_by(Group, ONTOLOGY) %>%
    # 选取 qvalue 最低的 10 个条目
    slice_min(order_by = qvalue, n = 10) %>% 
    ungroup() %>%
    # 为绘图准备 Description 因子
    mutate(Description = factor(Description, levels = unique(Description[order(qvalue)])))

# 绘制 GO 点图
if (nrow(plot_go_data) > 0) {
    p_go_dotplot <- ggplot(plot_go_data, aes(x = Count, y = Description, color = p.adjust, size = Count)) +
        geom_point() +
        # FIX: 将 ONT 更改为正确的列名 ONTOLOGY
        facet_grid(ONTOLOGY ~ Group, scales = "free") +
        scale_color_viridis_c(direction = -1) +
        theme_bw() +
        labs(
            title = "Top 10 GO Enrichment Analysis (K2-C1 vs K2-C2)",
            x = "Gene Count",
            y = "GO Term",
            color = "Adjusted P-value"
        ) +
        theme(strip.text.y = element_text(angle = 0))

    ggsave("07_Neutrophil_GO_Dotplot.png", p_go_dotplot, width = 12, height = 10)
    message("GO 富集点图已保存: 07_Neutrophil_GO_Dotplot.png")
} else {
    warning("GO 富集结果为空，跳过 GO 点图绘制。")
}

message("--- 07 脚本执行完毕 ---")