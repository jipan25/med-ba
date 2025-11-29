# R Script: geo/data_prep_ba.R
# 用于胆道闭锁 (Biliary Atresia) GEO 数据集 (GSE122340) 的差异表达分析
# 使用 GEOquery 和 limma 包，提供比 Python 更稳定的数据加载能力。

# --- 1. 安装和加载必要的包 ---
# 注意: 如果遇到 "compilation failed" 错误 (例如，缺少 x86_64-conda-linux-gnu-cc)，
# 请在您的 Linux 或 Conda 环境中安装编译器，例如:
# conda install -c conda-forge gcc gfortran

# 检查并安装 Bioconductor 包管理器
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# 检查并安装所需的 GEOquery 和 limma 包
required_packages <- c("Biobase", "GEOquery", "limma")

# 统一检查和安装所有包，以便 BiocManager 更好地处理依赖关系
packages_to_install <- required_packages[!sapply(required_packages, require, character.only = TRUE)]

# if (length(packages_to_install) > 0) {
if ( 0 > 1 )
    cat(paste0(" -> 正在安装缺失的包: ", paste(packages_to_install, collapse = ", "), "\n"))
    # 使用 BiocManager 一次性安装所有缺失的包及其依赖
    BiocManager::install(packages_to_install, update = FALSE, dependencies = TRUE)
}

# 确保所有包都已加载
sapply(required_packages, library, character.only = TRUE)

# --- 2. 配置参数 ---
GEO_ID <- "GSE122340"
OUTPUT_DIR <- "geo_analysis_results"
GROUP_KEY_SEARCH <- "disease state" # 用于在表型数据中搜索分组信息的关键词
CONTROL_TERM <- "normal"
CASE_TERM <- "Biliary atresia"
LFC_THRESHOLD <- 1.0  # Log2 Fold Change 阈值
P_VALUE_THRESHOLD <- 0.05 # 经过FDR校正后的P值 (adj.P.Val) 阈值

# 创建输出目录
if (!dir.exists(OUTPUT_DIR)) {
    dir.create(OUTPUT_DIR)
}

# --- 3. 主要分析函数 ---
run_dea <- function(geo_id, output_dir, group_key, control_term, case_term, lfc_thresh, p_thresh) {
    
    cat(paste0("\n--- 准备开始胆道闭锁生信复现之旅：GEO 数据的差异表达分析 (R/limma) ---\n"))
    cat(paste0("--- 尝试加载 GEO 数据集: ", geo_id, " ---\n"))

    # 1. 下载和加载 ExpressionSet
    tryCatch({
        # 使用 GSEMatrix=TRUE 直接下载 Series Matrix 文件，提高稳定性
        gset <- getGEO(geo_id, GSEMatrix = TRUE, Annot = TRUE, destdir=file.path(output_dir, geo_id))
        
        if (length(gset) > 1) {
            # 针对只有一个平台的情况
            idx <- grep("GPL", names(gset))
            gset <- gset[[idx]]
        } else {
            gset <- gset[[1]]
        }
    }, error = function(e) {
        cat(paste0("[ERROR] 加载 GEO 数据集失败。原始错误: ", e$message, "\n"))
        return(NULL)
    })
    
    if (is.null(gset)) return(NULL)

    cat(" -> 成功加载 ExpressionSet。\n")

    # 2. 提取表达矩阵和表型数据
    exprs_data <- exprs(gset)
    pData_data <- pData(gset)
    
    # 3. 识别分组信息
    group_col_name <- NULL
    
    # 尝试查找包含关键词的列
    for (col in colnames(pData_data)) {
        if (any(grepl(group_key, pData_data[[col]], ignore.case = TRUE))) {
            group_col_name <- col
            break
        }
    }
    
    if (is.null(group_col_name)) {
        cat("[ERROR] 无法在表型数据中找到包含关键词 'disease state' 的分组列。\n")
        return(NULL)
    }
    
    # 清理和创建因子
    groups <- pData_data[[group_col_name]]
    # 清理字符串，只保留关键词后的值 (e.g., "disease state: normal" -> "normal")
    groups <- gsub(paste0(".*", group_key, ":[[:space:]]*"), "", groups, ignore.case = TRUE)
    
    # 将分组转换为因子，确保对照组是第一个水平 (Reference)
    groups_factor <- factor(groups, levels = c(control_term, case_term))
    
    if (length(levels(groups_factor)) < 2 || any(is.na(groups_factor))) {
        cat("[ERROR] 分组信息不完整或只包含一组。请检查 CONTROL_TERM 和 CASE_TERM。\n")
        return(NULL)
    }
    
    cat(paste0(" -> 成功识别分组。控制组 ('", control_term, "') 样本数: ", sum(groups_factor == control_term), "\n"))
    cat(paste0(" -> 疾病组 ('", case_term, "') 样本数: ", sum(groups_factor == case_term), "\n"))

    # 4. 差异表达分析 (limma)
    cat("--- 开始进行差异表达分析 (limma) ---\n")
    
    # 创建设计矩阵
    design <- model.matrix(~0 + groups_factor)
    colnames(design) <- c(control_term, case_term)
    
    # 拟合线性模型
    fit <- lmFit(exprs_data, design)
    
    # 定义对比矩阵 (疾病组 vs. 控制组)
    contrast_matrix <- makeContrasts(
        Contrast = paste(case_term, control_term, sep = " - "),
        levels = design
    )
    
    fit2 <- contrasts.fit(fit, contrast_matrix)
    
    # 经验贝叶斯平滑
    fit2 <- eBayes(fit2)
    
    # 提取结果表 (按调整P值排序)
    # n=Inf 表示提取所有基因
    dea_results <- topTable(fit2, coef = "Contrast", number = Inf, adjust.method = "fdr")
    
    # 5. 清理和保存结果
    dea_results$ProbeID <- rownames(dea_results)
    
    # 过滤 DEGs
    degs <- subset(dea_results, 
                   adj.P.Val < p_thresh & 
                   abs(logFC) >= lfc_thresh)
                   
    degs$Regulation <- ifelse(degs$logFC > 0, "Up", "Down")

    # 完整结果路径
    full_path <- file.path(output_dir, paste0(geo_id, "_DEA_full_results.csv"))
    # DEGs 结果路径
    degs_path <- file.path(output_dir, paste0(geo_id, "_DEGs_LFC", lfc_thresh, "_FDR", p_thresh, ".csv"))

    # 保存完整结果
    write.csv(dea_results, full_path, row.names = FALSE)
    cat(paste0(" -> 完整DEA结果已保存到: ", full_path, "\n"))

    # 保存筛选后的DEGs
    write.csv(degs[c("ProbeID", "logFC", "adj.P.Val", "Regulation")], degs_path, row.names = FALSE)
    cat(paste0(" -> ", nrow(degs), " 个DEGs已保存到: ", degs_path, "\n"))
    cat(paste0(" -> 上调基因数: ", sum(degs$Regulation == 'Up'), "\n"))
    cat(paste0(" -> 下调基因数: ", sum(degs$Regulation == 'Down'), "\n"))
    
    # 打印DEGs前几行作为摘要
    if (nrow(degs) > 0) {
        cat("\n--- 差异表达基因 (Top 5) 摘要 ---\n")
        print(head(degs[c("ProbeID", "logFC", "adj.P.Val", "Regulation")], 5))
    } else {
        cat("\n--- 未发现满足阈值 (LFC>±1, FDR<0.05) 的差异表达基因 ---\n")
    }
    
    cat("\n--- 分析成功完成！请查看geo_analysis_results文件夹中的文件。 ---\n")
}

# 运行主函数
run_dea(GEO_ID, OUTPUT_DIR, GROUP_KEY_SEARCH, CONTROL_TERM, CASE_TERM, LFC_THRESHOLD, P_VALUE_THRESHOLD)