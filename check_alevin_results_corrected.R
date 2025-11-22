#!/usr/bin/env Rscript

library(Matrix)
library(ggplot2)

# 检查单个样本的函数
check_sample <- function(sample_name) {
  cat("=== 检查样本:", sample_name, "===\n")
  
  alevin_dir <- paste0("~/soft-dev/data/results/alevin/", sample_name)
  
  # 检查文件是否存在
  quants_file <- paste0(alevin_dir, "/alevin/quants_mat.gz")
  barcodes_file <- paste0(alevin_dir, "/alevin/quants_mat_cols.txt")
  genes_file <- paste0(alevin_dir, "/alevin/quants_mat_rows.txt")
  
  if (file.exists(quants_file) && file.exists(barcodes_file) && file.exists(genes_file)) {
    cat("✓ 找到所有必需文件\n")
    
    # 读取barcodes
    barcodes <- readLines(barcodes_file)
    cat("检测到的细胞barcode数量:", length(barcodes), "\n")
    
    # 读取基因名
    genes <- readLines(genes_file)
    cat("检测到的基因数量:", length(genes), "\n")
    
    # 读取量化矩阵
    # 注意：quants_mat.gz是稀疏矩阵格式
    counts <- readMM(gzfile(quants_file))
    
    # 设置行名和列名
    rownames(counts) <- genes
    colnames(counts) <- barcodes
    
    cat("\n量化矩阵维度:", dim(counts), "\n")
    cat("总UMI数:", sum(counts), "\n")
    
    # 计算质控指标
    umi_per_cell <- Matrix::colSums(counts)
    genes_per_cell <- Matrix::colSums(counts > 0)
    
    cat("\n每个细胞的UMI数统计:\n")
    print(summary(umi_per_cell))
    
    cat("\n每个细胞检测到的基因数统计:\n") 
    print(summary(genes_per_cell))
    
    # 检查线粒体基因比例
    mt_genes <- grep("^mt-", genes, ignore.case = TRUE, value = TRUE)
    if (length(mt_genes) > 0) {
      cat("找到", length(mt_genes), "个线粒体基因\n")
      mt_counts <- Matrix::colSums(counts[mt_genes, ])
      mt_percent <- mt_counts / umi_per_cell * 100
      cat("\n线粒体基因百分比统计:\n")
      print(summary(mt_percent))
      
      # 保存质控数据
      qc_data <- data.frame(
        cell_barcode = barcodes,
        nUMI = umi_per_cell,
        nGene = genes_per_cell,
        percent_mt = mt_percent
      )
    } else {
      cat("\n警告: 未找到线粒体基因\n")
      qc_data <- data.frame(
        cell_barcode = barcodes,
        nUMI = umi_per_cell,
        nGene = genes_per_cell
      )
    }
    
    return(qc_data)
    
  } else {
    cat("✗ 缺失必要文件\n")
    return(NULL)
  }
}

# 检查两个样本
qc_834 <- check_sample("CRR524834")
qc_835 <- check_sample("CRR524835")
