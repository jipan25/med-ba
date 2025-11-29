# 文件名: 2_01_download_datasets.R
# 项目: 2号研究 - BA 靶向细胞发现
# 功能: 自动从 GEO 数据库下载 GSE248744, GSE236230, GSE176189 的补充数据文件(scRNA-seq Raw Data)
# 环境: 此脚本设计用于 Docker 容器内执行，数据将保存在挂载的 /data 目录中。
# -----------------------------------------------------------------------------

# 1. 环境设置
# -----------------------------------------------------------------------------
# 增加下载超时时间（单细胞数据通常很大，默认为60秒可能不够）
options(timeout = 3600) 

# 检查并安装 GEOquery
if (!requireNamespace("GEOquery", quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos = "https://cloud.r-project.org")
    BiocManager::install("GEOquery", update = FALSE)
}
library(GEOquery)
library(tools) # 用于文件扩展名处理

# 2. 定义目标数据集和目录
# -----------------------------------------------------------------------------
# 定义数据集列表
# GSE248744: 胆管细胞异质性
# GSE236230: 肝细胞重编程
# GSE176189: 免疫微环境 & 空间转录组
target_gse_list <- c("GSE248744", "GSE236230", "GSE176189")

# !!! 关键路径修正 !!!
# 在 Docker 环境中，我们将数据保存到挂载的 /data 目录下。
# 对应本地路径为 /home/ma/data
data_root_dir <- "/data/BA_Study_2_Data"

# 创建主目录
if (!dir.exists(data_root_dir)) {
    dir.create(data_root_dir, recursive = TRUE) # recursive = TRUE 确保创建嵌套目录
    message(paste("已创建数据主目录 (位于 Docker 挂载点):", data_root_dir))
}

# 3. 批量下载函数
# -----------------------------------------------------------------------------
download_dataset <- function(gse_id, base_dir) {
    message(paste("\n========================================================"))
    message(paste("正在处理数据集:", gse_id))
    message(paste("========================================================"))
    
    # 检查目标文件夹是否存在（getGEOSuppFiles会自动创建）
    target_path <- file.path(base_dir, gse_id)
    if (dir.exists(target_path) && length(list.files(target_path)) > 0) {
        message(paste("注意: 目标目录", target_path, "已存在且包含文件。跳过下载。"))
        message("如果需要重新下载，请先手动删除该目录。")
        return()
    }
    
    tryCatch({
        # 1. 下载补充文件 (Raw Counts 通常在这里)
        message("正在尝试下载补充文件 (Supplementary Files)...")
        
        # getGEOSuppFiles 会自动在 base_dir 下创建一个以 gse_id 命名的文件夹
        files <- getGEOSuppFiles(gse_id, baseDir = base_dir, makeDirectory = TRUE)
        
        # 打印下载的文件列表
        message("下载成功! 文件列表:")
        print(rownames(files))
        
        # 2. 自动解压 (如果是 .tar 或 .tar.gz 包)
        dataset_dir <- file.path(base_dir, gse_id)
        downloaded_files <- list.files(dataset_dir, full.names = TRUE)
        
        for (file in downloaded_files) {
            # 检查是否是 tar 包
            if (grepl("\\.tar(\\.gz)?$", file, ignore.case = TRUE)) {
                message(paste("检测到压缩包，正在解压:", basename(file)))
                untar(file, exdir = dataset_dir)
                message("解压完成。")
            }
        }
        
    }, error = function(e) {
        message(paste("❌ 下载失败:", gse_id))
        message("错误信息:", e$message)
        message("提示: 如果网络连接不稳定，您可能需要手动从 NCBI GEO 网站下载，并放入目标目录中。")
    })
}

# 4. 执行下载
# -----------------------------------------------------------------------------
# 注意：我们传入 base_dir 的父目录，即 /data，让 GEOquery 在 /data 下创建 BA_Study_2_Data/GSE...
for (gse in target_gse_list) {
    # base_dir 应该是 /data
    download_dataset(gse, dirname(data_root_dir)) 
}

message("\n########################################################")
message("所有下载任务已尝试完成。")
message(paste("请检查本地路径 (Docker 容器外): /home/ma/data/BA_Study_2_Data"))
message(paste("请检查容器路径 (Docker 容器内):", data_root_dir))
message("下一步: 我们将检查文件格式并编写数据加载和 Seurat 对象创建脚本。")
message("########################################################")