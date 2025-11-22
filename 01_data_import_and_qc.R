# 文件名: 01_data_import_and_qc.R (修正版)
# 功能: 分步导入两个样本的 Salmon Alevin 数据，执行论文定义的严格 QC，然后合并。

library(Seurat)
library(tximport)
library(tidyverse)
library(EnsDb.Mmusculus.v79) 
library(Matrix)

# ==============================================================================
# 步骤 1: 配置路径和导入数据 (使用循环处理多样本)
# ==============================================================================

# --- A. 配置 Salmon Alevin 输出路径 ---
files <- c(
    "NC_5d" = "~/soft-dev/data/results/alevin_raw/CRR524834/alevin/quants_mat.gz",
    "RRV_5d" = "~/soft-dev/data/results/alevin_raw/CRR524835/alevin/quants_mat.gz"
)

if(!all(file.exists(files))) stop("错误：找不到部分 Salmon Alevin 结果文件。请检查路径。")

# --- B. 准备基因 ID 转换表 (Ensembl ID -> Gene Symbol) ---
message("--- 正在准备基因ID转换表 ---")
edb <- EnsDb.Mmusculus.v79
tx2gene <- transcripts(edb, return.type = "DataFrame")
tx2gene <- as.data.frame(tx2gene[,c("tx_id", "gene_id")])

# --- C. 循环导入和创建 Seurat 对象 ---
seurat_list <- list()

for (sample_name in names(files)) {
    file_path <- files[sample_name]
    message(paste("--- 正在导入样本:", sample_name, "---"))

    # 1. 导入数据 (tximport 一次只能处理一个 Alevin 文件)
    txi <- tximport(file_path, type = "alevin", tx2gene = tx2gene, ignoreTxVersion = TRUE)

    # 2. 创建 Seurat 对象
    sobj <- CreateSeuratObject(
        counts = txi$counts, 
        project = sample_name, 
        min.cells = 3
    )
    
    # 3. 添加样本信息到 metadata
    sobj$sample <- sample_name 

    message(paste(sample_name, "原始细胞数:", ncol(sobj)))

    # ==============================================================================
    # 步骤 2: 计算 QC 指标 (线粒体、血红蛋白)
    # ==============================================================================
    
    message("--- 正在计算 QC 指标 ---")

    # 1. 线粒体比例
    sobj[["percent.mt"]] <- PercentageFeatureSet(sobj, pattern = "^mt-")

    # 2. 血红蛋白比例 (论文要求的额外过滤)
    sobj[["percent.hb"]] <- PercentageFeatureSet(sobj, pattern = "^Hb[ab]-")

    # ==============================================================================
    # 步骤 3: 执行论文的严格质控 (QC)
    # ==============================================================================
    
    message("--- 正在执行论文的严格过滤标准 ---")
    # 论文标准: 200 < Genes < 5000; Mito < 5%; Hb < 5%
    sobj_filtered <- subset(sobj, subset = 
        nFeature_RNA > 200 & 
        nFeature_RNA < 5000 & 
        percent.mt < 5 & 
        percent.hb < 5
    )

    message(paste(sample_name, "过滤后细胞数:", ncol(sobj_filtered)))
    
    seurat_list[[sample_name]] <- sobj_filtered
}

# ==============================================================================
# 步骤 4: 合并过滤后的 Seurat 对象
# ==============================================================================

message("\n--- 正在合并过滤后的 Seurat 对象 ---")
# 使用 merge 函数将两个过滤后的对象合并成一个大对象
seurat_combined <- merge(seurat_list[[1]], y = seurat_list[2:length(seurat_list)], 
                         add.cell.ids = names(seurat_list), 
                         project = "Gr1_Cells_Combined")

message("==================================================================")
message(paste("最终合并后的细胞总数:", ncol(seurat_combined)))
message(paste("论文目标细胞总数:", 8351))
message("==================================================================")

# 保存合并且过滤后的对象，供下一步使用
saveRDS(seurat_combined, file = "seurat_combined_filtered.rds")
message("合并且过滤后的 Seurat 对象已保存为 seurat_combined_filtered.rds")

# 释放内存
rm(list = ls())