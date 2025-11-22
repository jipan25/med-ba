#!/usr/bin/env Rscript

library(Matrix)

# 读取单个样本的数据
read_alevin_data <- function(sample_name) {
  alevin_dir <- paste0("~/soft-dev/data/results/alevin/", sample_name)
  
  quants_file <- paste0(alevin_dir, "/alevin/quants_mat.gz")
  barcodes_file <- paste0(alevin_dir, "/alevin/quants_mat_cols.txt")
  genes_file <- paste0(alevin_dir, "/alevin/quants_mat_rows.txt")
  
  counts <- readMM(gzfile(quants_file))
  barcodes <- readLines(barcodes_file)
  genes <- readLines(genes_file)
  
  rownames(counts) <- genes
  colnames(counts) <- barcodes
  
  return(counts)
}

# 检查关键基因的表达
check_key_gene_expression <- function(counts, sample_name) {
  cat("\n=== 检查", sample_name, "中关键基因的表达 ===\n")
  
  # 论文相关的关键基因
  key_genes <- c(
    "Cd177",  # CD177+ cells - 论文重点
    "Elane",  # Neutrophil elastase - NETs相关
    "Mpo",    # Myeloperoxidase - NETs相关
    "Camp",   # Cathelicidin antimicrobial peptide - NETs相关
    "Lcn2",   # Lipocalin-2 - 中性粒细胞标志物
    "S100a8", # Calprotectin subunit - 中性粒细胞标志物
    "S100a9", # Calprotectin subunit - 中性粒细胞标志物
    "Cxcr2",  # Neutrophil chemokine receptor
    "Itgam"   # CD11b - 中性粒细胞标志物
  )
  
  # 在数据集中查找这些基因
  available_genes <- key_genes[key_genes %in% rownames(counts)]
  missing_genes <- key_genes[!key_genes %in% rownames(counts)]
  
  cat("找到的基因:", paste(available_genes, collapse = ", "), "\n")
  cat("缺失的基因:", paste(missing_genes, collapse = ", "), "\n")
  
  if (length(available_genes) > 0) {
    # 检查这些基因的表达
    gene_counts <- counts[available_genes, ]
    
    cat("\n关键基因的表达统计（UMI计数）:\n")
    for (gene in available_genes) {
      gene_expression <- gene_counts[gene, ]
      expressing_cells <- sum(gene_expression > 0)
      total_expression <- sum(gene_expression)
      
      cat(gene, ": ", 
          "在", expressing_cells, "个细胞中表达", 
          "（", round(expressing_cells/ncol(counts)*100, 2), "%），",
          "总UMI=", total_expression, "\n", sep = "")
    }
    
    # 特别关注CD177+细胞
    cd177_positive <- sum(counts["Cd177", ] > 0)
    if ("Cd177" %in% available_genes) {
      cat("\nCD177+ 细胞数量:", cd177_positive, 
          "（占所有细胞的", round(cd177_positive/ncol(counts)*100, 2), "%）\n")
      
      # 检查CD177+细胞的平均表达水平
      cd177_expression <- counts["Cd177", ]
      cd177_positive_cells <- cd177_expression > 0
      if (sum(cd177_positive_cells) > 0) {
        cat("CD177+细胞的平均CD177表达量:", 
            round(mean(cd177_expression[cd177_positive_cells]), 2), "UMI\n")
      }
    }
  }
  
  return(available_genes)
}

# 分析两个样本
cat("开始分析关键基因表达...\n")

counts_834 <- read_alevin_data("CRR524834")
genes_834 <- check_key_gene_expression(counts_834, "CRR524834")

counts_835 <- read_alevin_data("CRR524835") 
genes_835 <- check_key_gene_expression(counts_835, "CRR524835")
