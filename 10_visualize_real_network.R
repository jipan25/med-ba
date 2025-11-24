# 文件名: 10_visualize_real_network.R
# 功能: 读取 STRING 数据库的真实 PPI 数据，分析枢纽基因并可视化

library(tidyverse)
library(igraph)
library(ggraph)
library(AnnotationDbi)
library(EnsDb.Mmusculus.v79)

# ==============================================================================
# 1. 读取数据与自动诊断
# ==============================================================================
string_file <- "string_interactions.tsv"
if (!file.exists(string_file)) {
    stop("❌ 未找到 'string_interactions.tsv'。请先从 STRING 官网下载该文件并放入当前目录。")
}

message("--- 读取真实 PPI 网络数据 ---")
# 使用 read.delim 读取 TSV，更健壮
# quote="" 防止基因描述中的引号导致读取错误
ppi_data <- read.delim(string_file, header = TRUE, sep = "\t", quote = "")

message("--- 数据列名诊断 ---")
print(colnames(ppi_data))

# 智能查找 combined_score 列
# STRING 导出的列名可能是 combined_score, score, 或 score (combined_score) 等
score_col <- grep("score", colnames(ppi_data), value = TRUE, ignore.case = TRUE)

if (length(score_col) == 0) {
    stop("\n❌ 错误: 在文件中未找到 'combined_score' 或类似的得分列。\n",
         "请确认您下载的是 'TSV: tab separated values (lists only one-way edges)' 格式。\n",
         "当前文件的列名如上所示。")
}

# 如果找到多个（例如 experimental_score, combined_score），优先选 combined_score
if ("combined_score" %in% score_col) {
    target_score_col <- "combined_score"
} else {
    target_score_col <- score_col[1] # 否则取第一个带 score 的列
}
message(paste("✅ 锁定得分列:", target_score_col))

# 重命名该列为标准名称 combined_score，方便后续处理
colnames(ppi_data)[colnames(ppi_data) == target_score_col] <- "combined_score"

# 读取之前的 DEG 信息用于着色
message("--- 读取 DEG 信息 ---")
deg_info <- read.csv("06_Neutrophil_K2_DEG_results.csv")
# ID 清洗和 Symbol 映射
deg_info$clean_id <- gsub("\\.\\d+$", "", deg_info$gene)
deg_info$symbol <- mapIds(EnsDb.Mmusculus.v79, keys = deg_info$clean_id, column = "SYMBOL", keytype = "GENEID", multiVals = "first")

# ==============================================================================
# 2. 构建网络
# ==============================================================================
# 假设前两列是节点名称 (node1, node2)
# STRING 下载的数据通常前两列就是基因名
colnames(ppi_data)[1:2] <- c("node1", "node2")

# 过滤低置信度边 (combined_score > 0.4)
# 注意：STRING 的 combined_score 范围通常是 0-1，但也可能是 0-1000
# 我们先检查一下范围
max_score <- max(ppi_data$combined_score, na.rm = TRUE)
threshold <- 0.4
if (max_score > 1) {
    threshold <- 400 # 如果是 1000 分制，阈值设为 400
    message("检测到分数为 1000 分制，设置阈值为 400")
}

edges <- ppi_data %>%
    dplyr::filter(combined_score > threshold) %>%
    dplyr::select(node1, node2, combined_score)

if (nrow(edges) == 0) {
    stop("❌ 错误: 过滤后没有剩余的边。请检查阈值或数据文件内容。")
}

# 创建图对象
g <- graph_from_data_frame(d = edges, directed = FALSE)

# ==============================================================================
# 3. 计算枢纽基因 (Degree Centrality)
# ==============================================================================
V(g)$degree <- degree(g)

# 将 DEG 信息 (logFC, Group) 映射到网络节点
# 注意：STRING 导出的节点名就是基因 Symbol
# 使用 match 匹配
match_idx <- match(V(g)$name, deg_info$symbol)
V(g)$logFC <- deg_info$avg_log2FC[match_idx]
# 如果匹配不到（可能是STRING用了别名），默认为 NA
V(g)$Group <- ifelse(is.na(V(g)$logFC), "Interactor (Non-DEG)", 
                     ifelse(V(g)$logFC > 0, "K2-C1 Up (Antiviral)", "K2-C2 Up (Classical)"))

# 找出 Top Hubs
node_df <- data.frame(name = V(g)$name, degree = V(g)$degree, group = V(g)$Group) %>%
    arrange(desc(degree))

message("\n🏆 真实的 Top 10 枢纽基因 (Hub Genes):")
print(head(node_df, 10))

# 标记 Top 10 基因用于绘图
top_hubs <- head(node_df$name, 10)
V(g)$is_hub <- V(g)$name %in% top_hubs

# ==============================================================================
# 4. 可视化
# ==============================================================================
message("--- 正在绘制真实 PPI 网络 ---")

p_real <- ggraph(g, layout = 'fr') + 
    geom_edge_fan(aes(alpha = combined_score), color = "gray80", show.legend = FALSE) + 
    geom_node_point(aes(color = Group, size = degree), alpha = 0.9) +
    geom_node_text(
        aes(label = name), 
        data = . %>% dplyr::filter(name %in% top_hubs), 
        repel = TRUE, 
        fontface = "bold", 
        size = 4,
        color = "black"
    ) +
    scale_color_manual(
        values = c("K2-C1 Up (Antiviral)" = "#E41A1C", "K2-C2 Up (Classical)" = "#377EB8", "Interactor (Non-DEG)" = "gray60"),
        na.value = "gray60"
    ) +
    theme_graph() +
    labs(
        title = "Real PPI Network of Neutrophil Subsets",
        subtitle = "Data source: STRING Database | Top 10 Hubs Labeled"
    ) +
    theme(legend.position = "bottom")

ggsave("10_Real_Neutrophil_PPI_Network.png", p_real, width = 10, height = 9)
message("✅ 真实网络图已保存: 10_Real_Neutrophil_PPI_Network.png")