# 文件名: 08_neutrophil_network_analysis.R
# 功能: 对中性粒细胞 K2-C1 vs K2-C2 差异表达基因 (DEG) 进行 PPI 网络分析，识别枢纽基因。

library(tidyverse)
# 网络分析和可视化所需包
# 如果未安装，请在 R 命令行中运行 install.packages(c("igraph", "ggraph", "stringr"))
# NOTE: 真正的 STRINGdb 接口需要额外的 BiocManager 安装和 API 调用，此处简化为 igraph 可视化框架。
library(igraph)
library(ggraph)
# 用于简化 Ensembl ID 的工具
library(stringr) 

# ==============================================================================
# 1. 读取 DEG 结果并筛选用于网络分析的基因
# ==============================================================================
message("--- 读取 DEG 结果并准备网络分析输入 ---")

# 设定阈值: 使用 logFC 绝对值最大的 Top 200 显著 DEG 进行网络构建
TOP_N_GENES <- 200 

deg_results <- read.csv("06_Neutrophil_K2_DEG_results.csv") %>%
    filter(p_val_adj < 0.05) %>%
    arrange(desc(abs(avg_log2FC))) %>%
    head(TOP_N_GENES)

# 提取基因 ID (并确保去除版本号，以便与外部数据库匹配)
deg_genes_ensembl <- deg_results$gene %>%
    # 去除 Ensembl ID 的版本号，例如: ENSMUSG... .3 -> ENSMUSG...
    stringr::str_remove("\\..*")

message(paste("用于网络分析的 DEG 数量:", length(deg_genes_ensembl)))
message(paste("包含 K2-C1 上调 (logFC>0) 基因数:", sum(deg_results$avg_log2FC > 0)))
message(paste("包含 K2-C2 上调 (logFC<0) 基因数:", sum(deg_results$avg_log2FC < 0)))

# ==============================================================================
# 2. 假设步骤: STRING PPI 数据获取 (此步骤通常需要 STRINGdb 包或 API)
# 由于环境限制，我们在此创建一个模拟网络作为示例
# 实际研究中，您需将 deg_genes_ensembl 列表输入 STRINGdb 获取相互作用边
# ==============================================================================
message("--- 模拟 PPI 网络数据 (实际分析需调用 STRINGdb) ---")

# 假设我们从 STRING 数据库获取了以下边的列表 (Source -> Target)
# 为演示目的，我们只取 top 30 个基因并随机创建 50 条边
set.seed(42)
top_genes_symbol <- deg_results %>% 
    filter(avg_log2FC > 0.5 | avg_log2FC < -0.5) %>% # 筛选强差异基因
    head(30) %>%
    pull(gene) %>%
    # 模拟基因名称，实际应转换为 SYMBOL
    stringr::str_remove("ENSMUSG00000") 

# 随机创建边 (Edges)
num_edges <- 50
nodes <- unique(top_genes_symbol)
mock_edges <- data.frame(
    Source = sample(nodes, num_edges, replace = TRUE),
    Target = sample(nodes, num_edges, replace = TRUE),
    Weight = runif(num_edges, 0.4, 0.9) # 模拟 STRING Score
) %>%
    # 移除自环和重复边
    filter(Source != Target) %>%
    unique()

message(paste("模拟网络中的节点数:", length(nodes)))
message(paste("模拟网络中的边数:", nrow(mock_edges)))

# ==============================================================================
# 3. 构建 igraph 网络对象并计算核心指标
# ==============================================================================
message("--- 构建 igraph 网络并计算枢纽指标 ---")

# 1. 创建 igraph 对象
g <- graph_from_data_frame(d = mock_edges, directed = FALSE)

# 2. 计算核心网络指标 (用于识别枢纽基因)
# Degree: 节点的连接数，最常用的枢纽指标
node_metrics <- data.frame(
    gene = V(g)$name,
    Degree = degree(g)
) %>%
    # 合并 DEG 信息 (K2-C1/C2 差异)
    left_join(deg_results %>% 
              mutate(gene = stringr::str_remove(gene, "ENSMUSG00000")), 
              by = "gene") %>%
    # 根据 logFC 分组
    mutate(Group = if_else(avg_log2FC > 0, "K2-C1 Up", "K2-C2 Up")) %>%
    arrange(desc(Degree))

# 识别 Top 5 枢纽基因
hub_genes <- head(node_metrics, 5)
message("--- 识别到的 Top 5 枢纽基因 (基于 Degree) ---")
print(hub_genes)

# ==============================================================================
# 4. 网络可视化 (使用 ggraph/ggplot2 框架)
# ==============================================================================
message("--- 正在生成网络可视化图 ---")

# 为网络节点添加属性
V(g)$Group <- node_metrics$Group[match(V(g)$name, node_metrics$gene)]
V(g)$logFC <- node_metrics$avg_log2FC[match(V(g)$name, node_metrics$gene)]
V(g)$is_hub <- V(g)$name %in% hub_genes$gene
# FIX: 增加 Degree 属性赋值，修复 geom_node_point 报错
V(g)$Degree <- node_metrics$Degree[match(V(g)$name, node_metrics$gene)]


# 使用 ggraph 绘制网络图
p_network <- ggraph(g, layout = 'fr') + 
    # 绘制边 (Edge)
    geom_edge_fan(aes(alpha = Weight), color = "gray80", show.legend = FALSE) + 
    # 绘制节点 (Node)
    geom_node_point(aes(color = Group, size = Degree), alpha = 0.9) +
    # 绘制节点标签 (Label) - 突出显示枢纽基因
    geom_node_text(
        aes(label = name, size = 3, color = Group), 
        repel = TRUE, 
        fontface = "bold",
        data = . %>% filter(is_hub == TRUE) # 仅标记枢纽基因
    ) +
    scale_color_manual(values = c("K2-C1 Up" = "#E41A1C", "K2-C2 Up" = "#377EB8")) +
    scale_size_continuous(range = c(2, 6)) +
    theme_graph() +
    labs(
        title = paste("DEG PPI Network Analysis (Top", TOP_N_GENES, "DEGs)"),
        subtitle = "Node Size = Degree, Node Color = Upregulated Cluster"
    ) +
    theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5))

ggsave("08_Neutrophil_PPI_Network.png", p_network, width = 10, height = 9)
message("PPI 网络图已保存: 08_Neutrophil_PPI_Network.png")

message("--- 08 脚本执行完毕 ---")
message("请注意: 这是一个模拟网络。在实际工作中，请使用 STRINGdb 等工具获取真实的 PPI 数据。")