# 设置 CRAN 镜像源，提升下载速度
options(repos = c(CRAN = "https://cran.rstudio.com/"))

# --- 1. 安装 BiocManager ---
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager", quiet = TRUE)

# --- 2. 定义 Bioconductor 和 CRAN 包列表 ---

# Bioconductor 综合包列表
bioc_packages_comprehensive <- c(
    # 核心基础设施和数据结构
    "BiocGenerics", "IRanges", "S4Vectors", "GenomicRanges", "SummarizedExperiment",
    "Biobase", "AnnotationDbi", "DelayedArray", "XVector", "BiocVersion",

    # 转录组学/差异表达分析
    "DESeq2", "limma", "edgeR", "DEXSeq", "tximeta", "GenomicAlignments",

    # 单细胞 RNA-seq 分析 (scRNA-seq)
    "SingleCellExperiment", "scater", "scran", "DropletUtils", "scDblFinder",
    "celldex", "SingleR", "iSEE", "Rtsne", "umap",

    # 基因组学/表观遗传学 (ChIP-seq, ATAC-seq, 变异分析)
    "ChIPseeker", "DiffBind", "systemPipeR", "rtracklayer", "TxDb.Hsapiens.UCSC.hg38.knownGene",
    "BSgenome", "Gviz", "Rsamtools", "VariantAnnotation", "motifStack", "motifmatchr",

    # 注释和富集分析
    "org.Hs.eg.db", "GO.db", "KEGGREST", "clusterProfiler", "DOSE", "ReactomePA",
    "GSEABase", "GSVA",

    # 空间转录组学 (基础支持)
    "SpatialExperiment",

    # 高级可视化和报告
    "ComplexHeatmap", "pcaExplorer", "ggtree", "rmarkdown", "knitcitations", "bookdown"
)

# CRAN 拓展包列表 (添加 'future' 和 'furrr' 解决并行计算依赖问题)
cran_packages <- c(
    # 并行计算核心: 解决 'multisession' 错误
    "future",
    "furrr",

    # Tidyverse 全家桶 (数据清洗、处理、可视化)
    "tidyverse", "data.table", "reshape2", "magrittr", "glue",
    # 单细胞生态
    "Seurat", "SeuratObject", "sctransform",
    # 高级绘图
    "ggplot2", "pheatmap", "viridis", "RColorBrewer", "ggpubr", "cowplot",
    "patchwork", "plotly", "esquisse",
    # 实用工具
    "devtools", "usethis", "fs", "R.utils"
)

# --- 3. 执行安装 ---

message("--- 开始安装 Bioconductor 综合包集 (将跳过已存在的包) ---")
BiocManager::install(bioc_packages_comprehensive, update = FALSE, ask = FALSE, quiet = TRUE, force = TRUE)

message("--- 开始安装 CRAN 拓展包 ---")
install.packages(cran_packages, quiet = TRUE)

# --- 4. 安装 tinytex (用于 R Markdown 报告生成) ---

message("--- 检查并安装 tinytex (R Markdown 报告依赖) ---")
if (!requireNamespace("tinytex", quietly = TRUE)) {
    install.packages("tinytex", quiet = TRUE)
}
suppressMessages(library(tinytex))
if (!tinytex::is_tinytex()) {
    tinytex::install_tinytex(force = TRUE)
}

message("--- 所有包和依赖安装完成 ---")