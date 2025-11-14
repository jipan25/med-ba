#!/bin/bash

# 生物信息学工具安装脚本
# 为胆道闭锁分析安装必要的工具

echo "=== 安装生物信息学分析工具 ==="

# 检查conda是否安装
if ! command -v conda &> /dev/null; then
    echo "错误: 未找到conda，请先安装Miniconda或Anaconda"
    exit 1
fi

# 创建生物信息学环境
echo "创建生物信息学环境..."
conda create -n bioinfo python=3.9 -y

# 激活环境
source $(conda info --base)/etc/profile.d/conda.sh
conda activate bioinfo

# 安装核心生物信息学工具
echo "安装核心生物信息学工具..."
conda install -c bioconda -y \
    fastqc \
    multiqc \
    trim-galore \
    star \
    hisat2 \
    samtools \
    subread \
    htseq \
    deseq2 \
    r-base \
    bioconductor-edger

# 安装Python数据分析包
echo "安装Python数据分析包..."
pip install \
    pandas \
    numpy \
    scipy \
    matplotlib \
    seaborn \
    scikit-learn \
    jupyter

echo "=== 工具安装完成 ==="
echo ""
echo "使用方法:"
echo "1. 激活环境: conda activate bioinfo"
echo "2. 运行分析: python real_data_analysis_pipeline.py"
echo ""
echo "下一步:"
echo "1. 下载人类参考基因组"
echo "2. 配置基因组索引"
echo "3. 运行完整分析流程"