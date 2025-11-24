# 文件名: 09_prepare_string_input.R
# 功能: 提取 Top 200 DEG 的基因 Symbol，用于 STRING 数据库在线分析

library(tidyverse)
library(AnnotationDbi)
library(EnsDb.Mmusculus.v79)

# 1. 读取 DEG 结果
deg_results <- read.csv("06_Neutrophil_K2_DEG_results.csv")

# 2. 筛选 Top 200 显著基因
# FIX: 使用 dplyr::filter 避免与 ensembldb 冲突
top_deg <- deg_results %>%
    dplyr::filter(p_val_adj < 0.05) %>%
    arrange(desc(abs(avg_log2FC))) %>%
    head(200)

# 3. ID 转换 (Ensembl -> Symbol)
# 去除版本号
clean_ids <- gsub("\\.\\d+$", "", top_deg$gene)

# 映射 Symbol
symbols <- mapIds(EnsDb.Mmusculus.v79, 
                  keys = clean_ids, 
                  column = "SYMBOL", 
                  keytype = "GENEID", 
                  multiVals = "first")

# 4. 创建导出数据框
# FIX: 使用 dplyr::filter 和 dplyr::select
export_df <- top_deg %>%
    mutate(gene_symbol = symbols) %>%
    dplyr::filter(!is.na(gene_symbol)) %>% # 去除没找到名字的
    dplyr::select(gene_symbol, avg_log2FC)

# 5. 保存文件
write.table(export_df$gene_symbol, 
            "09_Gene_List_for_STRING.txt", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)

message("✅ 已生成基因列表文件: 09_Gene_List_for_STRING.txt")
message("请按照指导将此文件内容复制到 STRING 网站进行分析。")
message("\n--- 预览前 10 个基因 ---")
print(head(export_df$gene_symbol, 10))