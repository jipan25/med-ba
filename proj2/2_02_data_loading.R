# -------------------------------------------------------------------------
# 脚本名称: 2_04_data_loading.R
# 功能: 加载三个GEO数据集，进行基本QC，并保存为RDS文件。
# *** 修复 GSE176189 H5AD 转换逻辑 (强制重建统一 Assay) ***
# -------------------------------------------------------------------------

# 1. 加载必需的库
library(Seurat)
library(Matrix)
library(tidyverse)
library(edgeR) 
library(future)
library(furrr)
library(future.apply)
library(data.table) 
# --- 修正/新增的库 ---
library(SingleCellExperiment) 
library(zellkonverter)
library(SeuratObject) # V5 对象操作可能需要
# ---------------------

# 设置并行计算环境，以利用服务器多核优势
n_cores <- 8 
plan(multisession, workers = n_cores)
message(paste("已启用并行计算，使用", n_cores, "个核心。"))

# 设定数据路径和输出路径
DATA_DIR <- "/data/BA_Study_2_Data"
OUT_DIR <- "/ba/scRNA_Analysis"

# 确保输出目录存在
if (!dir.exists(OUT_DIR)) {
  dir.create(OUT_DIR, recursive = TRUE)
}

# =================================================================
# 辅助函数：加载并聚合 FeatureCounts 输出 (针对 GSE248744)
# =================================================================
load_count_files <- function(count_files) {
    
    message("正在并行读取并聚合 FeatureCounts 文件...")
    
    count_list <- future_lapply(count_files, function(f) {
        
        result <- tryCatch({
            
            data <- data.table::fread(
                f, 
                header = FALSE, 
                sep = "\t", 
                skip = 2, 
                data.table = FALSE, 
                select = c(1, 8) 
            )
            
            colnames(data) <- c("ID", "Count") 
            data$ID <- as.character(data$ID)
            data$Count <- as.integer(data$Count)
            data <- data[!startsWith(data$ID, "__"), ]
            
            if (any(duplicated(data$ID))) {
                data_agg <- aggregate(Count ~ ID, data = data, FUN = sum)
                data <- data_agg
            }
            
            rownames(data) <- data$ID
            data <- data[ , "Count", drop = FALSE] 
            sample_name <- gsub("_(RNA-seq_RAW264.7)?.count.txt.gz", "", basename(f))
            colnames(data) <- sample_name
            
            return(data)
            
        }, error = function(e) {
            warning(paste("读取文件失败 (GSE248744):", basename(f), "错误:", conditionMessage(e)))
            return(NULL) 
        })
        
        return(result)
        
    }, future.seed = TRUE)
    
    count_list <- Filter(Negate(is.null), count_list)
    
    if (length(count_list) == 0) {
      stop("所有 GSE248744 文件加载失败。")
    }
    
    count_matrix <- do.call(cbind, count_list)
    
    message("FeatureCounts 数据聚合和合并完成。")
    return(count_matrix)
}


# =================================================================
# 辅助函数：使用矩阵乘法处理稀疏矩阵重复基因
# =================================================================
aggregate_sparse_matrix <- function(counts_matrix) {
    gene_symbols <- rownames(counts_matrix)
    if (!any(duplicated(gene_symbols))) {
        return(counts_matrix)
    }
    
    unique_symbols <- unique(gene_symbols)
    
    # 1. 创建聚合矩阵 (Grouping Matrix)
    agg_matrix <- sparseMatrix(
        i = match(gene_symbols, unique_symbols),
        j = 1:length(gene_symbols),
        x = rep(1, length(gene_symbols)), 
        dims = c(length(unique_symbols), length(gene_symbols))
    )

    # 2. 矩阵乘法进行行聚合
    new_counts_matrix <- agg_matrix %*% counts_matrix
    
    # 3. 设置行名为唯一的基因符号
    rownames(new_counts_matrix) <- unique_symbols
    return(new_counts_matrix)
}


# -------------------------------------------------------------------------
# A. 处理 GSE248744 (Bulk RNA-seq)
# -------------------------------------------------------------------------
message("\n--- A. 正在处理 GSE248744 (Bulk RNA-seq) ---")
gse248744_dir <- file.path(DATA_DIR, "GSE248744")
count_files <- list.files(gse248744_dir, pattern = "*.count.txt.gz", full.names = TRUE)

if (length(count_files) > 0) {
  tryCatch({
    bulk_matrix <- load_count_files(count_files)
    bulk_matrix <- as.matrix(bulk_matrix)
    dge_list_248744 <- DGEList(counts = bulk_matrix)
    
    samples <- colnames(dge_list_248744$counts)
    pheno <- data.frame(SampleID = samples)
    
    pheno <- pheno %>%
      mutate(
        Condition = sub("GSM[0-9]+_RNA-seq_RAW264.7_(.*)_rep[0-9]+", "\\1", SampleID),
        Replicate = sub("GSM[0-9]+_RNA-seq_RAW264.7_.*_(rep[0-9]+)", "\\1", SampleID)
      )
    
    pheno$Condition <- gsub("_\\.", ".", pheno$Condition)
    pheno$Condition <- gsub("_", ".", pheno$Condition)
    pheno$Condition <- gsub("..siFACT", ".siFACT", pheno$Condition, fixed = TRUE)
    
    dge_list_248744$samples <- cbind(dge_list_248744$samples, pheno)
    
    message(paste("GSE248744 Bulk DGEList 对象创建成功。样本数:", ncol(dge_list_248744), "基因数:", nrow(dge_list_248744)))
    
    saveRDS(dge_list_248744, file.path(OUT_DIR, "GSE248744_DGEList.rds"))
    message(paste("GSE248744 DGEList 已保存到:", file.path(OUT_DIR, "GSE248744_DGEList.rds")))

  }, error = function(e) {
    warning(paste("GSE248744 处理失败，跳过。错误信息:", conditionMessage(e)))
  })
} else {
  warning("GSE248744 目录中未找到 *.count.txt.gz 文件，跳过处理。")
}


# -------------------------------------------------------------------------
# B. 处理 GSE236230 (scRNA-seq, 10X)
# -------------------------------------------------------------------------
message("\n--- B. 正在处理 GSE236230 (scRNA-seq, 10X) ---")
gse236230_dir <- file.path(DATA_DIR, "GSE236230")
samples_236230 <- c("GSM7520157_WT", "GSM7520158_TTP.KO")

seurat_list_236230 <- future_lapply(samples_236230, function(sample_name) {
  
  result <- tryCatch({
    files <- list(
      barcodes = file.path(gse236230_dir, paste0(sample_name, ".barcodes.tsv.gz")),
      features = file.path(gse236230_dir, paste0(sample_name, ".features.tsv.gz")),
      matrix = file.path(gse236230_dir, paste0(sample_name, ".matrix.mtx.gz"))
    )
    
    if (!all(file.exists(unlist(files)))) {
      message(paste("GSE236230 缺少 10X 必需文件:", sample_name, "跳过。"))
      return(NULL)
    }
    
    matrix_data <- readMM(files$matrix)
    features <- read.delim(files$features, header = FALSE, stringsAsFactors = FALSE)
    
    gene_symbols <- features$V2
    barcodes <- read.delim(files$barcodes, header = FALSE, stringsAsFactors = FALSE)
    colnames(matrix_data) <- barcodes$V1
    
    if (any(duplicated(gene_symbols))) {
        message(paste("GSE236230 -", sample_name, "发现", sum(duplicated(gene_symbols)), "个重复基因符号。正在进行聚合 (使用矩阵乘法)..."))
        counts_matrix <- aggregate_sparse_matrix(matrix_data)
    } else {
        rownames(matrix_data) <- gene_symbols
        counts_matrix <- matrix_data
    }
    
    seurat_obj <- CreateSeuratObject(
      counts = counts_matrix, 
      project = "GSE236230", 
      min.cells = 3, 
      min.features = 200
    )
    
    condition <- ifelse(grepl("WT", sample_name), "WT", "TTP.KO")
    seurat_obj$condition <- condition
    seurat_obj$sample_id <- sample_name
    
    seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^mt-")
    
    message(paste("GSE236230 -", sample_name, "对象创建成功。细胞数:", ncol(seurat_obj), "特征数:", nrow(seurat_obj)))
    return(seurat_obj)
  }, error = function(e) {
    warning(paste("GSE236230 样本", sample_name, "处理失败。错误:", conditionMessage(e)))
    return(NULL)
  })
  
  return(result)
  
}, future.seed = TRUE)

seurat_list_236230 <- Filter(Negate(is.null), seurat_list_236230)

if (length(seurat_list_236230) > 0) {
  tryCatch({
    sample_ids_for_merge <- purrr::map_chr(seurat_list_236230, ~unique(.x$sample_id))
    
    seurat_236230 <- merge(seurat_list_236230[[1]], 
                           y = seurat_list_236230[-1], 
                           add.cell.ids = sample_ids_for_merge, 
                           project = "GSE236230_Merged")

    saveRDS(seurat_236230, file.path(OUT_DIR, "GSE236230_Unfiltered.rds"))
    message(paste("GSE236230 合并对象已保存到:", file.path(OUT_DIR, "GSE236230_Unfiltered.rds")))
  }, error = function(e) {
    warning(paste("GSE236230 合并或保存失败。错误:", conditionMessage(e)))
  })
} else {
  warning("GSE236230 数据加载失败或无可用样本，跳过合并。")
}


# -------------------------------------------------------------------------
# C. 处理 GSE176189 (scRNA-seq + VDJ)
# -------------------------------------------------------------------------
message("\n--- C. 正在处理 GSE176189 (scRNA-seq + VDJ) ---")

gse176189_dir <- file.path(DATA_DIR, "GSE176189")

# 1. 确定 scRNA-seq 数据的目录
sc_10x_paths <- list.dirs(gse176189_dir, recursive = FALSE, full.names = TRUE) %>%
  keep(~grepl("filtered_feature_bc_matrix", basename(.x)))
  
if (length(sc_10x_paths) > 0) {
  
  message(paste("GSE176189 发现", length(sc_10x_paths), "个 10X 计数矩阵目录。"))
  
  # 2. 批量读取 scRNA-seq 数据
  sample_names_176189 <- basename(sc_10x_paths)
  
  seurat_list_176189 <- future_lapply(seq_along(sc_10x_paths), function(i) {
    path <- sc_10x_paths[i]
    sample_name <- sample_names_176189[i]
    
    result <- tryCatch({
        message(paste("-> 正在处理 GSE176189 样本:", sample_name))
        
        # 2a. 读取 10X 数据 (标准格式)
        data <- Read10X(data.dir = path)
        
        counts_matrix <- if (is.list(data)) { data$`Gene Expression` } else { data }
        
        if (is.null(counts_matrix) || prod(dim(counts_matrix)) == 0) {
          warning(paste("样本", sample_name, "的计数矩阵为空，跳过。"))
          return(NULL)
        }

        # 检查并聚合重复的基因符号
        if (any(duplicated(rownames(counts_matrix)))) {
            message(paste("GSE176189 -", sample_name, "发现重复基因符号。正在进行聚合 (使用矩阵乘法)..."))
            counts_matrix <- aggregate_sparse_matrix(counts_matrix)
        }

        seurat_obj <- CreateSeuratObject(
          counts = counts_matrix,
          project = "GSE176189", 
          min.cells = 3, 
          min.features = 200
        )
        
        # 2b. 添加元数据
        condition_id <- sub("([A-Z]+)_.*", "\\1", sample_name) 
        seurat_obj$condition_id <- condition_id
        seurat_obj$sample_group <- condition_id
        seurat_obj$sample_id <- sample_name
        
        # 2c. 整合 VDJ 数据 (TCR/BCR)
        vdj_files <- list.files(gse176189_dir, pattern = paste0("_", condition_id, "_(BCR|TCR)_filtered_contig_annotations.csv.gz"), full.names = TRUE, recursive = TRUE)
        
        if (length(vdj_files) > 0) {
          message(paste("GSE176189 -", condition_id, "找到 VDJ 文件，正在读取和整合..."))

          vdj_data_list <- lapply(vdj_files, function(f) {
            tryCatch({ data.table::fread(f, stringsAsFactors = FALSE) }, 
                     error = function(e) { 
                       warning(paste("读取 VDJ 文件失败:", basename(f), "错误:", conditionMessage(e)))
                       return(NULL)
                     })
          })
          
          vdj_data_list <- Filter(Negate(is.null), vdj_data_list)
          
          if (length(vdj_data_list) > 0) {
            full_vdj_dt <- data.table::rbindlist(vdj_data_list)
            
            vdj_meta <- as.data.frame(full_vdj_dt) %>%
              filter(is_cell == TRUE, high_confidence == TRUE) %>%
              dplyr::select(barcode, productive, chain, v_gene, j_gene, cdr3) %>%
              group_by(barcode) %>%
              summarise(
                VDJ_Chains = paste(unique(chain), collapse = ";"),
                VDJ_Productive = any(productive), 
                VDJ_CDR3s = paste(unique(cdr3[cdr3 != "None"]), collapse = ";"),
                .groups = 'drop'
              ) %>%
              column_to_rownames(var = "barcode")
            
            common_cells <- intersect(colnames(seurat_obj), rownames(vdj_meta))
            if (length(common_cells) > 0) {
                seurat_obj <- AddMetaData(seurat_obj, metadata = vdj_meta[common_cells, ])
            }
          }
        }
        
        message(paste("GSE176189 -", sample_name, "对象创建成功。细胞数:", ncol(seurat_obj)))
        return(seurat_obj)
        
    }, error = function(e) {
        warning(paste("GSE176189 样本", sample_name, "处理失败。错误:", conditionMessage(e)))
        return(NULL) 
    })
    
    return(result)
    
  }, future.seed = TRUE)

  seurat_list_176189 <- Filter(Negate(is.null), seurat_list_176189)

  # 合并 Seurat 对象
  if (length(seurat_list_176189) > 0) {
    tryCatch({
      sample_ids_for_merge <- purrr::map_chr(seurat_list_176189, ~unique(.x$sample_id))

      seurat_176189 <- merge(seurat_list_176189[[1]], 
                             y = seurat_list_176189[-1], 
                             add.cell.ids = sample_ids_for_merge, 
                             project = "GSE176189_Merged")

      # 3. QC 步骤
      seurat_176189[["percent.mt"]] <- PercentageFeatureSet(seurat_176189, pattern = "^MT-")

      # 4. 保存结果
      saveRDS(seurat_176189, file.path(OUT_DIR, "GSE176189_Unfiltered.rds"))
      message(paste("GSE176189 合并对象已保存到:", file.path(OUT_DIR, "GSE176189_Unfiltered.rds")))

    # 5. 转换为 AnnData (H5AD) 格式供 Python QC 脚本使用 -----------------------------
      h5ad_file_name <- "gse176189_raw_data.h5ad" # 2_03_gse176189_qc.py 要求的名称
      output_h5ad_file <- file.path(OUT_DIR, h5ad_file_name)

      message(paste("正在将 Seurat 对象转换为 AnnData (h5ad) 格式并保存为:", output_h5ad_file))

      # =================================================================
      # *** H5AD 转换优化和修复开始 (关键修复在此) ***
      # =================================================================
      
      message("--- H5AD 转换诊断信息 ---")
      message(paste("Seurat V5 对象维度 (Features x Cells):", paste(dim(seurat_176189), collapse = " x ")))
      message(paste("默认 Assay Layers (转换前):", paste(Layers(seurat_176189), collapse = ", ")))
      message(paste("细胞元数据行数:", nrow(seurat_176189@meta.data)))
      message("-------------------------")

      tryCatch({
          # 1. 强制获取 Seurat 合并后的统一计数矩阵
          # GetAssayData 应该返回对齐后的矩阵，即使 Layers() 列表混乱。
          # 默认情况下，V5 对象会尝试返回整个 Assay 的合并视图。
          full_counts_matrix <- GetAssayData(seurat_176189, assay = "RNA") 
          
          if (is.null(full_counts_matrix) || prod(dim(full_counts_matrix)) == 0) {
              stop("未能从合并的 Seurat 对象中提取到有效的统一计数矩阵。")
          }
          
          # 2. 创建一个干净的 RNA Assay 对象，只包含这个统一的 'counts' 层。
          # 这将完全覆盖掉所有混乱的 counts.X 子层。
          clean_rna_assay <- CreateAssay5(counts = full_counts_matrix)
          
          # 3. 将新的 Assay 赋值给 Seurat 对象
          seurat_176189[["RNA"]] <- clean_rna_assay
          
          message(paste("✅ RNA Assay 已清理并重建。当前 Layer:", paste(Layers(seurat_176189), collapse = ", ")))

          # 4. 创建 SingleCellExperiment (SCE) 对象 (使用清理后的 Seurat 对象)
          gse176189_sce_clean <- as.SingleCellExperiment(seurat_176189, assay = "RNA")
          
          # 5. 写入 H5AD 格式
          writeH5AD(gse176189_sce_clean, file = output_h5ad_file)
          
          message(paste("🎉 H5AD 文件已通过强制统一 Assay 结构成功创建:", output_h5ad_file))
          
      }, error = function(e) {
          # 如果连内部结构访问都失败，则报告最终错误
          stop(paste("GSE176189 H5AD 转换最终失败。错误:", conditionMessage(e)))
      })
      # ----------------------------------------------------------------------------------

    }, error = function(e) {
      warning(paste("GSE176189 合并或保存失败。错误:", conditionMessage(e)))
    })
  } else {
    warning("GSE176189 数据加载失败或无可用样本，跳过合并。")
  }
} else {
  warning("GSE176189 目录中未找到任何包含 'filtered_feature_bc_matrix' 的子目录，跳过处理。")
}

message("\n所有数据集加载和初始处理完成。")

# 报告所有生成的RDS文件
message(paste("\n处理后的数据已保存到以下文件：", 
              list.files(OUT_DIR, pattern = "*.rds", full.names = TRUE), 
              collapse = "\n"))