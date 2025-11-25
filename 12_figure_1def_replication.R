# 文件名: 12_figure_1def_replication.R
# 功能: 修正 Figure 1F：复现论文中显示样本来源叠加的 UMAP 图
# -----------------------------------------------------------------------------

library(ggplot2)
library(dplyr)
# 确保 ggpubr 包已安装
if (!requireNamespace("ggpubr", quietly = TRUE)) {
    stop("需要 'ggpubr' 包来组合图形 (ggarrange)。请在 R 控制台运行: install.packages('ggpubr')")
}
library(ggpubr) # 用于组合图形

# ==============================================================================
# 1. 读取数据
# ==============================================================================
message("--- 正在读取 real_umap_metadata.csv ---")
umap_data <- read.csv("real_umap_metadata.csv")

# 确保 cluster 列是因子，并设置颜色
umap_data$cluster <- factor(umap_data$cluster, levels = c("K2-C1", "K2-C2"))
# 确保 sample_group 列是因子
umap_data$sample_group <- as.factor(umap_data$sample_group)

# 定义颜色
cluster_colors <- c("K2-C1" = "#E41A1C", "K2-C2" = "#377EB8")
# 定义样本颜色 (对应论文中的 Saline 和 RRV 组)
# 假设 NC_5d 对应 Saline (红色), RRV_5d 对应 RRV (蓝色)
sample_colors <- c("NC_5d" = "#E41A1C", "RRV_5d" = "#377EB8") 

# 为了更好地重现 Figure 1F，我们只需要关注 K2-C1 和 K2-C2 亚群的细胞。
# -----------------------------------------------------------------------------


# ==============================================================================
# 2. 绘制 Figure 1D: UMAP 亚群分布图 (未更改)
# ==============================================================================

message("--- 正在绘制 Figure 1D: UMAP 亚群分布 ---")
p_umap_cluster <- ggplot(umap_data, aes(x = UMAP_1, y = UMAP_2, color = cluster)) +
    geom_point(size = 0.5, alpha = 0.8) +
    scale_color_manual(values = cluster_colors, name = "Subcluster") +
    theme_minimal() +
    theme(
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "right",
        plot.title = element_text(hjust = 0.5)
    ) +
    labs(title = "Neutrophil Subclusters (K=2)")


# ==============================================================================
# 3. 绘制 Figure 1E: Cd177 表达量 FeaturePlot (未更改)
# ==============================================================================

message("--- 正在绘制 Figure 1E: Cd177 表达量 ---")

max_exp <- max(umap_data$CD177_expression)
min_exp <- min(umap_data$CD177_expression)

umap_data_sorted_exp <- umap_data %>% 
  arrange(CD177_expression)

p_umap_cd177 <- ggplot(umap_data_sorted_exp, aes(x = UMAP_1, y = UMAP_2, color = CD177_expression)) +
    geom_point(size = 0.5, alpha = 0.8) +
    scale_color_gradientn(
        colors = c("lightgrey", "yellow", "red"), # 灰-黄-红 梯度
        limits = c(min_exp, max_exp),
        name = "Cd177\nExpr."
    ) +
    theme_minimal() +
    theme(
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "right",
        plot.title = element_text(hjust = 0.5)
    ) +
    labs(title = "Cd177 Expression")


# ==============================================================================
# 4. 绘制 Figure 1F (修正): UMAP 样本来源叠加图 (类似论文中的图)
# ==============================================================================

message("--- 正在绘制 Figure 1F (修正): UMAP 样本来源叠加图 ---")

# 按照 sample_group 排序，将其中一组（比如 RRV_5d/蓝色）放在上层，以更好地进行可视化。
umap_data_sorted_sample <- umap_data %>% 
  arrange(sample_group)

p_umap_sample_overlay <- ggplot(umap_data_sorted_sample, aes(x = UMAP_1, y = UMAP_2, color = sample_group)) +
    geom_point(size = 0.5, alpha = 0.8) +
    scale_color_manual(values = sample_colors, name = "Sample Group") +
    theme_minimal() +
    theme(
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "right",
        plot.title = element_text(hjust = 0.5)
    ) +
    labs(title = "Sample Origin Overlay")

# ==============================================================================
# 5. 组合图形
# ==============================================================================

message("--- 正在组合 Figure 1 D, E, F ---")

# D (亚群) 和 E (表达量)
p_combined_de <- ggarrange(p_umap_cluster, p_umap_cd177, 
                           ncol = 2, 
                           labels = c("D", "E"),
                           common.legend = FALSE)

# 将组合图 (D+E) 和 F (样本来源) 垂直组合
p_final <- ggarrange(p_combined_de, p_umap_sample_overlay, 
                     ncol = 1, 
                     heights = c(1, 1), 
                     labels = c("", "F")) 

# 保存最终组合图
ggsave("Figure_1DEF_Combined_Replication_RevisedF.png", p_final, width = 10, height = 10)
message("🎉 最终组合图已保存为 Figure_1DEF_Combined_Replication_RevisedF.png (F图已修正为 UMAP 叠加图)")
message("现在 Figure 1F 应该与论文图相似，展示了亚群中不同样本的分布。")