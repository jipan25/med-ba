# 文件名: 12_figure_1def_replication.R
# 功能: 使用 05 脚本导出的 CSV 文件复现论文 Figure 1 D, E, F
# 修正: 将 Figure 1F 修改为 UMAP 样本来源叠加图，以匹配论文原图
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
# 这个文件由 05_neutrophil_subanalysis.R 脚本生成
data_file_path <- "real_umap_metadata.csv"

if (!file.exists(data_file_path)) {
    stop("错误: 找不到 'real_umap_metadata.csv' 文件。请确保已成功运行 05_neutrophil_subanalysis.R。")
}

umap_data <- read.csv(data_file_path)

# 确保 cluster 列是因子，并设置颜色
umap_data$cluster <- factor(umap_data$cluster, levels = c("K2-C1", "K2-C2"))
# 确保 sample_group 列是因子
umap_data$sample_group <- as.factor(umap_data$sample_group)

# 定义颜色
# 图 D: 亚群颜色
cluster_colors <- c("K2-C1" = "#1f78b4", "K2-C2" = "#ff7f00") # 蓝色/橙色 (参考论文配色)
# 图 F: 样本颜色 (Saline=红色, RRV=蓝色)
sample_colors <- c("NC_5d" = "#E41A1C", "RRV_5d" = "#377EB8") # 请根据您的实际样本名调整 Key

# 定义绘图主题
theme_custom <- function() {
  theme_minimal(base_size = 14) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black"),
      panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
      legend.title = element_text(face = "bold"),
      axis.text = element_blank(), # 移除轴刻度
      axis.ticks = element_blank(),
      axis.title = element_blank() # 移除轴标题 (如 UMAP_1)
    )
}

# ==============================================================================
# 2. 绘制 Figure 1D: UMAP 亚群分布图
# ==============================================================================

message("--- 正在绘制 Figure 1D: UMAP 亚群分布 ---")
p_umap_cluster <- ggplot(umap_data, aes(x = UMAP_1, y = UMAP_2, color = cluster)) +
    geom_point(size = 0.5, alpha = 0.8) +
    scale_color_manual(values = cluster_colors, name = "Subcluster") +
    labs(title = "D") +
    theme_custom() +
    theme(plot.title = element_text(size = 20, face = "bold", hjust = 0)) # 左对齐大标题

# ==============================================================================
# 3. 绘制 Figure 1E: Cd177 表达量 FeaturePlot
# ==============================================================================

message("--- 正在绘制 Figure 1E: Cd177 表达量 ---")

# 排序以防遮挡
umap_data_sorted_exp <- umap_data %>% arrange(CD177_expression)

p_umap_cd177 <- ggplot(umap_data_sorted_exp, aes(x = UMAP_1, y = UMAP_2, color = CD177_expression)) +
    geom_point(size = 0.5, alpha = 0.8) +
    scale_color_gradientn(
        colors = c("lightgrey", "magenta"), # 仿照论文图E的紫色
        name = "CD177"
    ) +
    labs(title = "E") +
    theme_custom() +
    theme(plot.title = element_text(size = 20, face = "bold", hjust = 0))

# ==============================================================================
# 4. 绘制 Figure 1F (修正): UMAP 样本来源叠加图
# ==============================================================================

message("--- 正在绘制 Figure 1F (修正): UMAP 样本来源叠加图 ---")

# 随机打乱顺序，避免一种颜色完全覆盖另一种
set.seed(123)
umap_data_shuffled <- umap_data[sample(nrow(umap_data)), ]

# 检查样本名是否匹配颜色定义的 Key
print("样本组名称:")
print(unique(umap_data_shuffled$sample_group))
# 如果您的样本名不是 NC_5d/RRV_5d，请在这里动态调整 sample_colors 的 names

p_umap_sample <- ggplot(umap_data_shuffled, aes(x = UMAP_1, y = UMAP_2, color = sample_group)) +
    geom_point(size = 0.5, alpha = 0.8) +
    scale_color_manual(values = sample_colors, name = "Group") +
    labs(title = "F") +
    theme_custom() +
    theme(plot.title = element_text(size = 20, face = "bold", hjust = 0))

# ==============================================================================
# 5. 组合图形
# ==============================================================================

message("--- 正在组合 Figure 1 D, E, F ---")

# 垂直排列 D, E, F
p_final <- ggarrange(p_umap_cluster, p_umap_cd177, p_umap_sample,
                     ncol = 1, 
                     nrow = 3,
                     align = "v")

# 保存最终组合图
ggsave("Figure_1DEF_Combined_Replication_RevisedF.png", p_final, width = 6, height = 15, dpi = 300)
message("🎉 最终组合图已保存为 Figure_1DEF_Combined_Replication_RevisedF.png")