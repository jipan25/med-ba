#!/usr/bin/env Rscript

library(Matrix)
library(ggplot2)
library(patchwork)

# 读取数据函数
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
  
  return(list(counts = counts, barcodes = barcodes, genes = genes))
}

# 读取两个样本
data_834 <- read_alevin_data("CRR524834")
data_835 <- read_alevin_data("CRR524835")

# 创建质控图
create_qc_plots <- function(data, sample_name) {
  counts <- data$counts
  
  # 计算QC指标
  qc_data <- data.frame(
    cell_barcode = colnames(counts),
    nUMI = Matrix::colSums(counts),
    nGene = Matrix::colSums(counts > 0)
  )
  
  # 计算线粒体百分比
  mt_genes <- grep("^mt-", rownames(counts), ignore.case = TRUE, value = TRUE)
  if (length(mt_genes) > 0) {
    mt_counts <- Matrix::colSums(counts[mt_genes, ])
    qc_data$percent_mt <- mt_counts / qc_data$nUMI * 100
  }
  
  # 创建图表
  p1 <- ggplot(qc_data, aes(x = nUMI)) +
    geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7) +
    labs(title = paste(sample_name, "- UMI分布"), 
         x = "每个细胞的UMI数", y = "细胞数量") +
    theme_minimal()
  
  p2 <- ggplot(qc_data, aes(x = nGene)) +
    geom_histogram(bins = 50, fill = "darkorange", alpha = 0.7) +
    labs(title = paste(sample_name, "- 基因数分布"),
         x = "每个细胞检测到的基因数", y = "细胞数量") +
    theme_minimal()
  
  if ("percent_mt" %in% colnames(qc_data)) {
    p3 <- ggplot(qc_data, aes(x = percent_mt)) +
      geom_histogram(bins = 50, fill = "forestgreen", alpha = 0.7) +
      labs(title = paste(sample_name, "- 线粒体基因百分比"),
           x = "线粒体基因百分比", y = "细胞数量") +
      theme_minimal()
    
    p4 <- ggplot(qc_data, aes(x = nUMI, y = nGene, color = percent_mt)) +
      geom_point(alpha = 0.5, size = 0.5) +
      scale_color_gradient(low = "blue", high = "red", name = "MT%") +
      labs(title = paste(sample_name, "- UMI vs 基因数"),
           x = "UMI数", y = "基因数") +
      theme_minimal()
    
    return((p1 + p2) / (p3 + p4))
  } else {
    return(p1 + p2)
  }
}

# 生成质控图
p834 <- create_qc_plots(data_834, "CRR524834")
p835 <- create_qc_plots(data_835, "CRR524835")

# 保存图表
ggsave("CRR524834_qc_plots.png", p834, width = 12, height = 8, dpi = 300)
ggsave("CRR524835_qc_plots.png", p835, width = 12, height = 8, dpi = 300)

cat("已生成质量控制图表\n")
cat("样本CRR524834: ", ncol(data_834$counts), "个细胞, ", nrow(data_834$counts), "个基因\n")
cat("样本CRR524835: ", ncol(data_835$counts), "个细胞, ", nrow(data_835$counts), "个基因\n")
