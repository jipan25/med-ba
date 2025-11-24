# 文件名: 11_figure_1ghi_replication.R
# 功能: 复现论文 Figure 1G, 1H, 1I
# 核心逻辑: 在 CD177+ 亚群内部，比较 RRV(感染) vs NC(对照)

library(Seurat)
library(tidyverse)
library(ggplot2)
library(clusterProfiler)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(EnsDb.Mmusculus.v79)
library(cowplot) # 用于组合图形

# ==============================================================================
# 1. 数据准备：提取 CD177+ 亚群 (K2-C1)
# ==============================================================================
message("--- 读取中性粒细胞子集数据 ---")
neutrophils <- readRDS("neutrophils_subclustered.rds")

# 关键修正：设置默认 Assay 为 RNA
DefaultAssay(neutrophils) <- "RNA"

# 确认当前的亚群名称
# 使用模糊匹配查找 Cd177，避免版本号问题
all_genes <- rownames(neutrophils)
cd177_id <- grep("ENSMUSG00000029657", all_genes, value = TRUE)[1]
if(is.na(cd177_id)) cd177_id <- grep("^Cd177$", all_genes, value = TRUE)[1]

if(is.na(cd177_id)) stop("无法找到 Cd177 基因，无法确定目标亚群。")

# 计算平均表达量以确定哪个簇是 CD177+
# 使用 LayerData 代替 AverageExpression 以避免警告，且更直接
expr_mat <- GetAssayData(neutrophils, assay = "RNA", layer = "data")
gene_expr <- expr_mat[cd177_id, ]
avg_expr_by_cluster <- tapply(gene_expr, neutrophils$sub_cluster, mean)
target_cluster <- names(avg_expr_by_cluster)[which.max(avg_expr_by_cluster)]

message(paste("检测到 CD177 高表达亚群为:", target_cluster))
message("--- 正在提取该亚群进行疾病组对比 ---")

# 提取 CD177+ 细胞
cd177_pos_cells <- subset(neutrophils, subset = sub_cluster == target_cluster)

# 设置对比身份为 'sample' (RRV_5d vs NC_5d)
Idents(cd177_pos_cells) <- "sample"

message("样本分布:")
print(table(Idents(cd177_pos_cells)))

# ==============================================================================
# 2. Figure 1G: 寻找 RRV 组上调最显著的基因
# ==============================================================================
message("--- 执行差异分析: RRV_5d vs NC_5d (in CD177+ cells) ---")

# 检查是否有足够的细胞进行比较
if(min(table(Idents(cd177_pos_cells))) < 3) {
  stop("某一组细胞数量太少，无法进行差异分析。")
}

deg_rrv <- FindMarkers(
    cd177_pos_cells,
    ident.1 = "RRV_5d",
    ident.2 = "NC_5d",
    verbose = FALSE,
    logfc.threshold = 0.25,
    min.pct = 0.1
)

# 添加基因 Symbol 列
deg_rrv$gene_id <- rownames(deg_rrv)
deg_rrv$clean_id <- gsub("\\.\\d+$", "", deg_rrv$gene_id)
deg_rrv$symbol <- mapIds(EnsDb.Mmusculus.v79, keys = deg_rrv$clean_id, column = "SYMBOL", keytype = "GENEID", multiVals = "first")
# 如果没有映射到，保留原ID
deg_rrv$symbol[is.na(deg_rrv$symbol)] <- deg_rrv$gene_id[is.na(deg_rrv$symbol)]

# 筛选上调基因
up_genes <- deg_rrv %>% 
    dplyr::filter(avg_log2FC > 0, p_val_adj < 0.05) %>%
    dplyr::arrange(desc(avg_log2FC))

message(paste("RRV 组显著上调基因数量:", nrow(up_genes)))

# 取 Top 8 (复现 Fig 1G)
top8_genes <- head(up_genes, 8)
message("--- Top 8 RRV Upregulated Genes (Fig 1G candidates) ---")
print(top8_genes[, c("symbol", "avg_log2FC", "p_val_adj")])

# --- 绘图逻辑 (纯 ggplot2 绘制) ---
message("--- 正在绘制 Fig 1G 小提琴图 (纯 ggplot2) ---")

# 1. 准备绘图特征列表
features_to_plot <- top8_genes$gene_id
# 使用 Symbol 作为列表名称，用于后续重命名
names(features_to_plot) <- top8_genes$symbol

# 2. 从 Seurat 对象中提取数据
vars_to_fetch <- c("sample", features_to_plot)
plot_data <- Seurat::FetchData(
    cd177_pos_cells, 
    vars = vars_to_fetch, 
    layer = "data" # 确保使用 data 层 (log-normalized expression)
)

# 3. 重命名列名，方便绘图 (使用 Symbol)
colnames(plot_data)[-1] <- names(features_to_plot)

# 4. 转换数据为长格式 (用于 facet_wrap 绘图)
plot_data_long <- plot_data %>%
  tidyr::pivot_longer(
    cols = -sample, # 排除 sample 列
    names_to = "Gene",
    values_to = "Expression"
  )

# 5. 使用 ggplot2 绘制小提琴图，并分面
p_1g <- ggplot(plot_data_long, aes(x = sample, y = Expression, fill = sample)) +
  # 绘制小提琴图，并设置透明度
  geom_violin(scale = "width", alpha = 0.8) +
  # 添加箱线图 (可选，但能显示中位数和四分位数)
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white", color = "black") + 
  # 添加原始数据点，使用抖动和低透明度
  geom_point(position = position_jitter(width = 0.15), size = 0.2, alpha = 0.4) + 
  # 分面，每行 4 个图，y 轴独立缩放
  facet_wrap(~ Gene, scales = "free_y", ncol = 4) +
  # 颜色设置
  scale_fill_manual(values = c("NC_5d" = "#A8DADC", "RRV_5d" = "#E63946")) +
  labs(title = "Top 8 RRV Upregulated Genes (CD177+ Neutrophils)", y = "Log Normalized Expression") +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    strip.text = element_text(face = "italic", size = 10, colour = "black"), # 确保基因名清晰
    panel.grid.major.x = element_blank(), # 移除垂直网格线
    panel.grid.minor = element_blank()
  )

# 6. 保存图形
ggsave("11_Fig1G_Top8_Upregulated.png", p_1g, width = 12, height = 6)
message("✅ Fig 1G 已保存。")

# 检查论文提到的特定基因
target_literatures <- c("Isg15", "Ifit1", "Ifi27l2a")
found_literatures <- up_genes %>% dplyr::filter(symbol %in% target_literatures)
if(nrow(found_literatures) > 0) {
    message("✅ 成功在差异列表中找到论文提到的关键基因:")
    print(found_literatures[, c("symbol", "avg_log2FC")])
} else {
    warning("⚠️ 未在前列差异基因中找到 Isg15, Ifit1 或 Ifi27l2a。")
}

# ==============================================================================
# 3. Figure 1H: 富集分析 (Metascape 替代方案)
# ==============================================================================
message("--- 执行功能富集分析 (复现 Fig 1H) ---")

# 转换 ID 为 Entrez
entrez_ids <- mapIds(org.Mm.eg.db, keys = up_genes$clean_id, column = "ENTREZID", keytype = "ENSEMBL", multiVals = "first")
entrez_ids <- na.omit(entrez_ids)

if(length(entrez_ids) > 0) {
    # GO BP 分析
    ego <- enrichGO(
        gene = entrez_ids,
        OrgDb = org.Mm.eg.db,
        ont = "BP",
        pAdjustMethod = "BH",
        qvalueCutoff = 0.05,
        readable = TRUE
    )
    
    if(!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
        # 筛选论文提到的关键词
        keywords <- "virus|interferon|antigen|response to external stimulus"
        selected_go <- as.data.frame(ego) %>% dplyr::filter(grepl(keywords, Description, ignore.case = TRUE))
        
        if (nrow(selected_go) > 0) {
            message("✅ 发现与论文描述一致的通路:")
            print(head(selected_go$Description, 5))
        }
        
        # 优化 dotplot 显示，只显示 Q 值最低的
        p_1h <- dotplot(ego, showCategory = 15) + ggtitle("Enrichment: RRV vs NC (CD177+)")
        ggsave("11_Fig1H_Enrichment.png", p_1h, width = 10, height = 8)
        message("✅ Fig 1H 已保存。")
    } else {
        warning("未发现显著富集的 GO 通路。")
    }
} else {
    warning("没有有效的 Entrez ID 用于富集分析。")
}

# ==============================================================================
# 4. Figure 1I: 功能热图 (Functional Differences)
# ==============================================================================
message("--- 执行功能热图分析 (复现 Fig 1I) ---")

# 使用 Top 30 上调基因作为默认或回退
features_to_heatmap <- head(up_genes$gene_id, 30)
heatmap_title <- "Top 30 Upregulated Genes in RRV"

# 绘制热图
# 关键：先对子集进行 ScaleData，指定要画图的基因
cd177_pos_cells <- ScaleData(cd177_pos_cells, features = features_to_heatmap, assay = "RNA", verbose = FALSE)

p_1i <- DoHeatmap(cd177_pos_cells, features = features_to_heatmap, group.by = "sample", label = FALSE) +
        scale_fill_gradientn(colors = c("blue", "white", "red")) +
        ggtitle(heatmap_title) + 
        theme(axis.text.y = element_text(size = 6)) # 缩小基因名以适应

ggsave("11_Fig1I_Functional_Heatmap.png", p_1i, width = 8, height = 10)
message("✅ Fig 1I 已保存。")

# ==============================================================================
# 5. 脚本收尾: 强制关闭图形设备 (解决零长度变量错误)
# ==============================================================================
if (!is.null(dev.list())) {
  dev.off() 
}

message("--- 11 脚本执行完毕 ---")