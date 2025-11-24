# 文件名: 10.2_visualize_real_network.R
# 功能: 模拟单细胞数据，并复现论文中 UMAP 聚类图 (D)、基因表达图 (E) 和分组分布图 (F)。
# **已修复分组操作中的尺寸不匹配错误**

# ------------------------------ 1. 加载必要的库 ------------------------------
# 确保您已安装这些库：install.packages(c("ggplot2", "dplyr", "patchwork"))
library(ggplot2)
library(dplyr)
library(patchwork) # 用于将多个图表组合在一起

# ------------------------------ 2. 模拟数据生成 ------------------------------
# 模拟 1000 个 Gr-1+ 细胞的数据点
n_cells <- 1000

# K-means 聚类 (K=2) 比例: C1 (~76.3%) 和 C2 (~23.7%)
n_c1 <- round(n_cells * 0.763) # 763
n_c2 <- n_cells - n_c1          # 237

set.seed(42) # 设置随机种子以保证结果可重复

# --- 2.1 模拟 K2-C1 (主体) 数据 ---
# 模拟 C1 (主体) 的 UMAP 坐标
df_c1 <- data.frame(
  UMAP_1 = rnorm(n_c1, mean = 0, sd = 2),
  UMAP_2 = rnorm(n_c1, mean = 0, sd = 2)
) %>%
  mutate(
    cluster = factor("K2-C1"),
    # 模拟 C1 的样本分组：假设 C1 中 Saline 占多数 (70% Saline, 30% RRV)
    sample_group = factor(
      sample(c("Saline", "RRV"), size = n_c1, replace = TRUE, prob = c(0.7, 0.3))
    ),
    # 模拟 CD177 表达：C1 簇低表达
    CD177_expression = rpois(n_c1, lambda = 0.5)
  )

# --- 2.2 模拟 K2-C2 (尾巴) 数据 ---
# 模拟 C2 (分离的尾巴) 的 UMAP 坐标
df_c2 <- data.frame(
  UMAP_1 = rnorm(n_c2, mean = 5, sd = 0.8), # 向右上方偏移
  UMAP_2 = rnorm(n_c2, mean = -4, sd = 0.8) # 向右下方偏移
) %>%
  mutate(
    cluster = factor("K2-C2"),
    # 模拟 C2 的样本分组：强制 RRV 富集 (25% Saline, 75% RRV)
    sample_group = factor(
      sample(c("Saline", "RRV"), size = n_c2, replace = TRUE, prob = c(0.25, 0.75))
    ),
    # 模拟 CD177 表达：C2 簇高表达
    CD177_expression = rpois(n_c2, lambda = 4) + 1
  )


# --- 2.3 合并数据并创建 CD177+ 标记 ---
umap_data <- bind_rows(df_c1, df_c2) %>%
  mutate(
    CD177_is_positive = ifelse(CD177_expression > 1, "CD177+", "Other")
  )

# ------------------------------ 3. 绘图复现 (图 D) ------------------------------
# 图 D: UMAP 聚类图，按 K-means 结果上色
plot_d <- ggplot(umap_data, aes(x = UMAP_1, y = UMAP_2, color = cluster)) +
  geom_point(size = 0.5) +
  scale_color_manual(
    values = c("K2-C1" = "#1f78b4", "K2-C2" = "#ff7f00"), # 对应原图的蓝色和橙色
    name = NULL
  ) +
  # 添加聚类百分比标签 (使用模拟的实际比例)
  annotate("text", x = -3.5, y = 3.5, label = paste0(round(n_c1/n_cells*100, 1), "%"), size = 6, color = "#1f78b4") +
  annotate("text", x = 6.5, y = -1.5, label = paste0(round(n_c2/n_cells*100, 1), "%"), size = 6, color = "#ff7f00") +
  labs(
    title = "D",
    x = "UMAP_1",
    y = "UMAP_2",
    caption = "Gr-1+"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 28, face = "bold", hjust = 0),
    legend.position = "right",
    legend.text = element_text(size = 14),
    axis.title = element_blank(), # 隐藏坐标轴标题
    axis.text = element_blank(),   # 隐藏坐标轴刻度文本
    axis.ticks = element_blank(),  # 隐藏坐标轴刻度线
    panel.grid = element_blank(),
    plot.caption = element_text(hjust = 0.5, size = 18)
  ) +
  coord_fixed() # 保持坐标轴比例一致

# ------------------------------ 4. 绘图复现 (图 E) ------------------------------
# 图 E: UMAP 基因表达图，按 CD177 表达水平上色
plot_e <- ggplot(umap_data, aes(x = UMAP_1, y = UMAP_2, color = CD177_expression)) +
  geom_point(size = 0.5) +
  # 使用紫色渐变来突出表达
  scale_color_gradientn(
    colors = c("gray90", "gray50", "#9400D3"), # 从灰色到紫色的渐变
    values = c(0, 0.5, 1),
    limits = c(0, max(umap_data$CD177_expression)),
    name = "CD177 Expression"
  ) +
  # 突出显示 CD177+ 细胞区域 (原图只显示了 CD177+ 细胞，这里用表达量模拟)
  geom_point(data = filter(umap_data, CD177_is_positive == "CD177+"), 
             color = "#9400D3", size = 0.5) +
  # 标签表示 CD177+ 细胞的比例
  annotate("text", x = -3.5, y = 3.5, label = "70.9%", size = 6, color = "#9400D3") +
  labs(
    title = "E",
    x = "UMAP_1",
    y = "UMAP_2",
    caption = "Gr-1+"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 28, face = "bold", hjust = 0),
    legend.position = "none", # 隐藏图例以模拟原图样式
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    plot.caption = element_text(hjust = 0.5, size = 18)
  ) +
  coord_fixed()

# ------------------------------ 5. 绘图复现 (图 F) ------------------------------
# 图 F: UMAP 样本分组分布图，按 RRV 和 Saline 样本上色
umap_data_f <- umap_data %>%
  # 重新排列绘图顺序，让 RRV 组点在上面 (模拟 RRV 组富集)
  arrange(sample_group)

plot_f <- ggplot(umap_data_f, aes(x = UMAP_1, y = UMAP_2, color = sample_group)) +
  geom_point(size = 0.5) +
  scale_color_manual(
    values = c("Saline" = "#e31a1c", "RRV" = "#1f78b4"), # 对应原图的红色和蓝色
    name = NULL
  ) +
  # 添加组别分布百分比标签 (使用原图 F 的数字来标注，模拟其位置)
  annotate("text", x = -3.5, y = 3.5, label = "66.7%", size = 6, color = "#1f78b4") + # RRV 组在主体部分的比例
  annotate("text", x = 6.5, y = -1.5, label = "74.5%", size = 6, color = "#e31a1c") + # Saline 组在尾巴部分的比例
  labs(
    title = "F",
    x = "UMAP_1",
    y = "UMAP_2",
    caption = "CD177+"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 28, face = "bold", hjust = 0),
    legend.position = "right",
    legend.text = element_text(size = 14),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    plot.caption = element_text(hjust = 0.5, size = 18)
  ) +
  coord_fixed()

# ------------------------------ 6. 组合并保存图片 ------------------------------
# 使用 patchwork 将 D, E, F 三个图组合成一个大图
combined_plot <- (plot_d | plot_e | plot_f) + 
  plot_layout(widths = c(1, 1, 1))

# 打印组合图
print(combined_plot)

ggsave("figure_1_d_e_f_combined.png", plot = combined_plot, width = 15, height = 5, units = "in", dpi = 300)