# 文件名: 06_neutrophil_deg_analysis.R
# 功能: 对中性粒细胞 K2-C1 和 K2-C2 亚群进行差异表达基因 (DEG) 分析

library(Seurat)
library(tidyverse)
library(ggplot2)

# 确保安装了用于火山图的包
if (!requireNamespace("ggrepel", quietly = TRUE)) {
    # 假设用户环境允许安装
    # install.packages("ggrepel", repos = "http://cran.us.r-project.org") 
    message("请确保已安装 ggrepel 包以显示基因标签。")
}
library(ggrepel)

# ==============================================================================
# 1. 读取数据并设置身份
# ==============================================================================
message("--- 读取中性粒细胞子集对象 ---")
# This file is saved by the 05 script
neutrophils <- readRDS("neutrophils_subclustered.rds")

# 设置默认 Assay 和 身份 (Idents)
DefaultAssay(neutrophils) <- "RNA"
Idents(neutrophils) <- "sub_cluster"

message("--- 检查亚群身份 ---")
print(table(Idents(neutrophils)))

# ==============================================================================
# 2. 差异表达基因 (DEG) 分析
# ==============================================================================
message("--- 执行 DEG 分析: K2-C1 (Majority) vs K2-C2 (Minority) ---")

# 使用 FindMarkers 进行比较:
# ident.1: K2-C1 (Majority) (高表达的基因将是 K2-C1 的标志基因)
# ident.2: K2-C2 (Minority)
# min.pct = 0.25, logfc.threshold = 0.25 是 Seurat 默认的保守阈值
deg_results <- FindMarkers(
    object = neutrophils,
    ident.1 = "K2-C1 (Majority)",
    ident.2 = "K2-C2 (Minority)",
    test.use = "wilcox", # 使用 Wilcoxon Rank Sum test，Seurat 默认
    min.pct = 0.25,
    logfc.threshold = 0.25,
    only.pos = FALSE # 包含上调和下调基因
)

# 清理结果并添加基因名称
deg_results <- deg_results %>%
    rownames_to_column(var = "gene") %>%
    # 筛选显著性基因
    filter(p_val_adj < 0.05) %>%
    # 根据 log2FC 绝对值排序
    arrange(desc(abs(avg_log2FC)))

message(paste("发现显著差异表达基因 (p_adj < 0.05, |logFC| > 0.25) 共:", nrow(deg_results)))

# 保存 DEG 结果
write.csv(deg_results, file = "06_Neutrophil_K2_DEG_results.csv", row.names = FALSE)
message("DEG 结果已保存至: 06_Neutrophil_K2_DEG_results.csv")

# ==============================================================================
# 3. 可视化 - 火山图 (Volcano Plot)
# ==============================================================================
message("--- 生成火山图 ---")

# 准备火山图数据
plot_data <- deg_results %>%
    mutate(
        # -log10(p_val_adj)
        logP = -log10(p_val_adj),
        # 定义基因状态
        Status = case_when(
            avg_log2FC > 0.25 & p_val_adj < 0.05 ~ "Up (K2-C1)",
            avg_log2FC < -0.25 & p_val_adj < 0.05 ~ "Down (K2-C2)",
            TRUE ~ "N.S."
        )
    )

# 标记前 15 个最显著的基因 (高logFC)
top_genes_to_label <- plot_data %>%
    filter(Status != "N.S.") %>%
    arrange(desc(abs(avg_log2FC))) %>%
    head(15)

p_volcano <- ggplot(plot_data, aes(x = avg_log2FC, y = logP, color = Status, label = gene)) +
    geom_point(alpha = 0.8, size = 1.5) +
    # 显著性阈值线
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray60") +
    geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed", color = "gray60") +
    scale_color_manual(values = c("Up (K2-C1)" = "#E41A1C", "Down (K2-C2)" = "#377EB8", "N.S." = "gray")) +
    # 添加基因标签，避免重叠
    geom_text_repel(
        data = top_genes_to_label,
        aes(label = gene),
        size = 3,
        box.padding = unit(0.5, "lines"),
        point.padding = unit(0.5, "lines"),
        max.overlaps = 20 # 允许更多重叠以显示关键基因
    ) +
    theme_minimal() +
    labs(
        title = "Volcano Plot: K2-C1 vs K2-C2 Neutrophils",
        x = "log2(Fold Change) [K2-C1 / K2-C2]",
        y = "-log10(Adjusted P-value)"
    ) +
    theme(legend.position = "bottom")

ggsave("06_Neutrophil_K2_Volcano_Plot.png", p_volcano, width = 8, height = 7)
message("火山图已保存: 06_Neutrophil_K2_Volcano_Plot.png")

# ==============================================================================
# 4. 可视化 - 热图 (Heatmap)
# ==============================================================================
message("--- 生成 DEG 热图 ---")

# 选择 Top 10 上调和 Top 10 下调基因
top_deg_up <- deg_results %>%
    filter(avg_log2FC > 0) %>%
    head(10) %>%
    pull(gene)

top_deg_down <- deg_results %>%
    filter(avg_log2FC < 0) %>%
    head(10) %>%
    pull(gene)

top_deg_genes <- unique(c(top_deg_up, top_deg_down))

# ------------------------------------------------------------------------------
# 修复步骤: 对用于绘制热图的基因子集进行 ScaleData
# ------------------------------------------------------------------------------
if (length(top_deg_genes) > 1) {
    message(paste("正在对", length(top_deg_genes), "个 Top DEG 进行 ScaleData..."))
    # 在 RNA Assay 上对 selected features 进行缩放
    neutrophils <- ScaleData(
        object = neutrophils,
        features = top_deg_genes,
        assay = "RNA",
        verbose = FALSE
    )

    # 绘制热图
    p_heatmap <- DoHeatmap(
        object = neutrophils,
        features = top_deg_genes,
        group.by = "sub_cluster",
        label = TRUE,
        disp.min = -2, disp.max = 2 # 限制颜色范围
    ) +
    # 使用 Bioconda 提供的颜色方案
    scale_fill_gradientn(colors = c("blue", "white", "red")) +
    theme(
        axis.text.y = element_text(size = 8),
        plot.title = element_text(hjust = 0.5)
    ) +
    ggtitle("Top 20 DEGs between K2-C1 and K2-C2")
    
    ggsave("06_Neutrophil_K2_Heatmap.png", p_heatmap, width = 10, height = 8)
    message("热图已保存: 06_Neutrophil_K2_Heatmap.png")
} else {
    warning("差异表达基因数量不足，跳过热图绘制。")
}

message("--- 06 脚本执行完毕 ---")