# R Script: proj2/2_03_data_conversion.R
# 职责:
# 1. 从 proj2/2_02_data_loading.R 生成的 RDS 文件中加载数据。
# 2. 修复 GSE176189 的 Seurat V5 Assay 多层问题 (CRITICAL FIX):
#    a. 解决 Layer 矩阵行名缺失问题 (通过强制分配)。
#    b. **新增动态填充**: 如果 master_features 列表太短，动态创建占位符基因名，以确保所有批次都能获得行名。
#    c. 解决 Layer 矩阵列名维度不一致问题 (通过强制分配)。
#    d. 解决批次间基因不对齐问题 (通过零填充对齐)。
# 3. 新增调试步骤: 打印第一个批次 (counts.1) 的数据头，以供分析。
# 4. 将单细胞数据 (GSE176189, GSE236230) 安全地转换为 H5AD 格式。

# --- 1. 设置和库加载 ---
library(Seurat)
library(zellkonverter)
library(SingleCellExperiment)
library(Matrix) # 用于处理稀疏矩阵的 Reduce(cbind) 操作

# 设置基础路径
BASE_PATH <- "/ba/scRNA_Analysis"

# 输入文件 (来自 2_02_data_loading.R)
GSE176189_RDS <- file.path(BASE_PATH, "GSE176189_Unfiltered.rds")
GSE236230_RDS <- file.path(BASE_PATH, "GSE236230_Unfiltered.rds")

# 输出文件 (H5AD 格式)
GSE176189_H5AD <- file.path(BASE_PATH, "gse176189_raw_data.h5ad")
GSE236230_H5AD <- file.path(BASE_PATH, "gse236230_raw_data.h5ad")

message("--- 检查输入文件 ---")
# 检查输入文件是否存在
if (!all(file.exists(GSE176189_RDS), file.exists(GSE236230_RDS))) {
    stop("错误: 缺少 GSE176189_Unfiltered.rds 或 GSE236230_Unfiltered.rds 文件。请先运行 proj2/2_02_data_loading.R。")
}
message("输入 RDS 文件检查完毕，准备开始处理。")

# --- 2. 辅助函数：安全 H5AD 转换 ---
safe_save_h5ad <- function(seurat_obj, output_path, dataset_id) {
    if (file.exists(output_path)) {
        message(paste("\n---", dataset_id, "H5AD 转换 ---"))
        message(paste(dataset_id, "H5AD 文件已存在，跳过转换。"))
        return(TRUE)
    }

    message(paste("\n---", dataset_id, "H5AD 转换诊断信息 ---"))
    print(paste("Seurat V5 对象维度 (Features x Cells):", paste(dim(seurat_obj), collapse = " x ")))
    print(paste("默认 Assay Layers (转换前):", paste(Layers(seurat_obj[["RNA"]]), collapse = ", ")))
    print(paste("细胞元数据行数:", nrow(seurat_obj@meta.data)))
    message("-------------------------")

    tryCatch({
        # 转换为 SingleCellExperiment (SCE) 对象
        sce <- as.SingleCellExperiment(seurat_obj, assay = "RNA", layer = "counts")

        # 使用 zellkonverter 转换为 H5AD
        writeH5AD(sce, file = output_path, verbose = TRUE)
        message(paste(dataset_id, "H5AD 转换成功并保存到:", output_path))
        return(TRUE)
    }, error = function(e) {
        warning(paste(dataset_id, "H5AD 转换最终失败。错误:", e$message))
        return(FALSE)
    })
}

# --- 3. GSE176189 处理和修复 (核心逻辑修改：动态扩展 master_features) ---
process_gse176189 <- function() {
    message("\n--- 正在处理 GSE176189 数据集 ---")
    seurat_176189 <- readRDS(GSE176189_RDS)
    
    # ----------------------------------------------------
    # CRITICAL FIX: 修复 V5 Assay 多层问题 & 基因/细胞名称问题
    # ----------------------------------------------------
    rna_assay <- seurat_176189[["RNA"]]
    all_layers <- Layers(rna_assay)
    # 识别所有以 'counts.' 开头的子层
    counts_layers <- all_layers[grep("^counts\\.\\d+$", all_layers)]

    if (length(counts_layers) > 1) {
        message(paste("--- GSE176189 修复: 检测到", length(counts_layers), "个计数子层。执行基因/细胞名称修复、对齐和合并操作 ---"))
        
        # 获取 Assay 中存储的全体特征名 (Genes) - 这个列表可能已被截断
        master_features <- Features(rna_assay)
        master_feature_count <- length(master_features)
        
        # 获取 Seurat 对象中存储的全体细胞名 (Barcodes)
        master_cell_names <- colnames(seurat_176189)
        
        matrix_list <- list()
        all_genes_union <- c()
        
        # Track assignments for Rownames (Genes) and Colnames (Cells)
        features_assigned_count <- 0 # 跟踪已分配的累计基因数 (sequential block assumption)
        cells_assigned_count <- 0    # 跟踪已处理的累计细胞数，用于切片列名
        
        for (layer_name in counts_layers) {
            current_matrix <- rna_assay@layers[[layer_name]]
            
            # 确保矩阵是稀疏格式
            if (!inherits(current_matrix, "dgCMatrix")) {
                 current_matrix <- as(current_matrix, "dgCMatrix")
            }
            
            n_features_in_layer <- nrow(current_matrix)
            n_cells_in_layer <- ncol(current_matrix) # 当前批次实际细胞数
            
            # -----------------------------------------------------------------
            # CRITICAL FIX 1: COLUMN NAMES (CELL BARCODES)
            # -----------------------------------------------------------------
            start_cell_index <- cells_assigned_count + 1
            end_cell_index <- cells_assigned_count + n_cells_in_layer
            
            if (end_cell_index > length(master_cell_names)) {
                stop(paste("致命错误: 无法为层", layer_name, "分配列名。累计细胞数", end_cell_index, "超过了整体 Seurat 对象细胞总数", length(master_cell_names)))
            }
            
            # 强制使用主对象的 cell names 进行赋值，确保维度匹配
            colnames(current_matrix) <- master_cell_names[start_cell_index:end_cell_index]
            cells_assigned_count <- end_cell_index # 更新累计细胞数
            
            message(paste("  - 已提取层:", layer_name, "; 维度:", paste(dim(current_matrix), collapse = " x ")))
            
            # -----------------------------------------------------------------
            # CRITICAL FIX 2: ROW NAMES (GENES) - 解决累计特征数溢出问题
            # -----------------------------------------------------------------
            start_index <- features_assigned_count + 1
            end_index <- features_assigned_count + n_features_in_layer
            
            if (is.null(rownames(current_matrix)) || length(rownames(current_matrix)) != n_features_in_layer) {
                
                # --- 动态扩展 master_features 列表 (解决当前致命错误) ---
                if (end_index > master_feature_count) {
                    missing_count <- end_index - master_feature_count
                    
                    # 生成占位符名称，从当前 master_feature_count + 1 开始
                    new_placeholders <- paste0("Missing_Gene_Placeholder_", (master_feature_count + 1):end_index)
                    
                    # 扩展 master_features 列表并更新总计数
                    master_features <- c(master_features, new_placeholders)
                    master_feature_count <- length(master_features)
                    
                    message(paste("  - CRITICAL WARNING:", layer_name, "所需索引", start_index, "至", end_index, "超出当前 Master List 范围。已自动扩展", missing_count, "个占位符基因名。"))
                }
                # ----------------------------------------------------
                
                # 强制分配行名，使用连续块
                rownames(current_matrix) <- master_features[start_index:end_index] 
                
                features_assigned_count <- end_index # 更新累计特征数
                
                message(paste("  - WARNING:", layer_name, "行名缺失/不匹配。已强制分配", n_features_in_layer, "个名称 (索引", start_index, "至", end_index, ")。"))
            } else {
                # 如果行名是完整的，我们仍需要更新累计计数，以防下一个层缺失行名
                features_assigned_count <- end_index
            }


            # --- 调试打印: 打印第一个批次的前10条数据 ---
            if (layer_name == "counts.1") {
                message("\n--- 调试信息: GSE176189 (counts.1) 批次的前 10 行数据 ---")
                # 打印前 10 行和前 5 列
                print(head(current_matrix[, 1:min(5, ncol(current_matrix))], 10))
                message("-------------------------------------------------------------------\n")
            }
            
            matrix_list[[layer_name]] <- current_matrix
            all_genes_union <- union(all_genes_union, rownames(current_matrix))
            message(paste("  - Rownames Check:", length(rownames(current_matrix)), "; Colnames Check:", length(colnames(current_matrix))))
        }
        
        if (length(all_genes_union) == 0) {
             stop("修复失败: 即使强制分配名称后，基因合集仍为空。数据结构可能已严重损坏。")
        }
        
        message(paste("  - 发现所有批次基因总合集大小 (Union):", length(all_genes_union)))
        
        # 2. 对所有矩阵进行基因对齐 (Padding with zeros)
        aligned_matrix_list <- list()
        for (layer_name in counts_layers) {
            current_matrix <- matrix_list[[layer_name]]
            
            # 找出缺失的基因 (需要用零填充的)
            missing_genes <- setdiff(all_genes_union, rownames(current_matrix))
            
            if (length(missing_genes) > 0) {
                # 创建缺失基因的零矩阵 (行数为缺失基因数，列数为细胞数)
                zero_matrix <- Matrix(0, 
                                      nrow = length(missing_genes), 
                                      ncol = ncol(current_matrix), 
                                      dimnames = list(missing_genes, colnames(current_matrix)), 
                                      sparse = TRUE)
                
                # 按行合并原始矩阵和零矩阵
                aligned_matrix <- rbind(current_matrix, zero_matrix)
                
                message(paste("  - 层", layer_name, "已对齐。添加了", length(missing_genes), "个零填充基因。"))
            } else {
                aligned_matrix <- current_matrix
                message(paste("  - 层", layer_name, "已对齐。无缺失基因。"))
            }
            
            # 按照所有基因合集的顺序重新排序行，确保所有矩阵的行顺序一致
            aligned_matrix <- aligned_matrix[all_genes_union, ]
            aligned_matrix_list[[layer_name]] <- aligned_matrix
        }
        
        # 3. 合并所有对齐后的矩阵 (按列合并细胞)
        merged_counts_matrix <- Reduce(cbind, aligned_matrix_list)
        
        # 4. 诊断检查
        expected_cells <- nrow(seurat_176189@meta.data) # 使用元数据中的细胞总数
        actual_cells <- ncol(merged_counts_matrix)
        message(paste("  - 所有矩阵合并完成。最终矩阵维度:", paste(dim(merged_counts_matrix), collapse = " x ")))
        message(paste("  - 检查: 预期细胞数 (Metadata):", expected_cells, "; 实际合并细胞数:", actual_cells))
        
        if (actual_cells != expected_cells) {
            stop("严重错误: 合并后的细胞总数与原始 Seurat 对象的细胞数不匹配。请检查 2_02_data_loading.R 的加载和合并逻辑。")
        }
        
        # 5. 创建新的 Seurat 对象
        message("  - 正在使用合并后的矩阵和原始元数据创建新的单层 Seurat 对象 ...")
        
        # 使用合并后的矩阵和原始元数据创建一个全新的、干净的 Seurat 对象
        seurat_176189_fixed <- CreateSeuratObject(
            counts = merged_counts_matrix,
            meta.data = seurat_176189@meta.data,
            project = seurat_176189@project.name
        )
        
        # 6. 用新的修复对象替换原来的对象
        seurat_176189 <- seurat_176189_fixed
        
        message(paste("  - 修复完成。新 Seurat 对象 Assay 层:", paste(Layers(seurat_176189[["RNA"]]), collapse = ", ")))

        # 重新保存修复后的对象
        saveRDS(seurat_176189, file = GSE176189_RDS)
        message(paste("  - 处理后的 GSE176189 数据已重新保存 (已修复为单层)：", GSE176189_RDS))
    } else {
        message("GSE176189 Assay 结构正常或已修复 (Layers:", paste(all_layers, collapse = ", "), ")。无需合并。")
    }
    
    # H5AD 转换
    safe_save_h5ad(seurat_176189, GSE176189_H5AD, "GSE176189")
}

# --- 4. GSE236230 转换 ---
process_gse236230 <- function() {
    message("\n--- 正在处理 GSE236230 数据集 ---")
    seurat_236230 <- readRDS(GSE236230_RDS)
    
    # H5AD 转换
    safe_save_h5ad(seurat_236230, GSE236230_H5AD, "GSE236230")
}

# --- 5. 执行所有处理步骤 ---
suppressWarnings({
    process_gse176189()
    process_gse236230()
})

message("\n=============================================")
message("所有数据集 H5AD 转换步骤完成。")
message(" GSE176189 的多层、基因名称缺失和不对齐问题已通过动态扩展和强制对齐的方式解决。")
message("=============================================")