#!/bin/bash
# install_conda_and_tools.sh

set -euo pipefail

echo "=== Conda及生物信息学工具安装脚本 ==="

# 检查是否已安装conda
if command -v conda &> /dev/null; then
    echo "✓ Conda 已安装"
else
    echo "安装Miniconda..."
    
    # 下载Miniconda
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda_installer.sh
    
    # 安装
    bash miniconda_installer.sh -b -p $HOME/miniconda3
    
    # 初始化conda
    $HOME/miniconda3/bin/conda init bash
    
    # 重新加载配置
    source ~/.bashrc
    
    echo "✓ Miniconda 安装完成"
fi

# 配置conda
echo "配置Conda..."
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict

# 创建质控专用环境
ENV_NAME="strict_qc_env"
echo "创建环境: $ENV_NAME"

if conda env list | grep -q "$ENV_NAME"; then
    echo "环境 $ENV_NAME 已存在，更新工具..."
    conda activate "$ENV_NAME"
else
    echo "创建新环境..."
    conda create -n "$ENV_NAME" -y python=3.8
    conda activate "$ENV_NAME"
fi

# 安装必要工具
echo "安装生物信息学工具..."

# 基础工具
conda install -y -c bioconda \
    trimmomatic=0.39 \
    fastqc=0.11.9 \
    multiqc=1.11 \
    fastp=0.23.2 \
    cutadapt=4.1 \
    seqtk=1.3

# 安装MultiQC（通过pip确保最新版）
pip install multiqc

# 验证安装
echo ""
echo "=== 验证安装 ==="

tools=("trimmomatic" "fastqc" "multiqc" "fastp" "cutadapt")
for tool in "${tools[@]}"; do
    if command -v "$tool" &> /dev/null; then
        version=$($tool --version 2>/dev/null | head -1 || echo "版本未知")
        echo "✓ $tool: $version"
    else
        echo "❌ $tool 未正确安装"
    fi
done

# 测试Trimmomatic适配器文件
echo ""
echo "=== 检查Trimmomatic适配器文件 ==="

ADAPTER_PATHS=(
    "$CONDA_PREFIX/share/trimmomatic/adapters/TruSeq3-PE-2.fa"
    "$(dirname $(which trimmomatic))/../share/trimmomatic/adapters/TruSeq3-PE-2.fa"
)

for path in "${ADAPTER_PATHS[@]}"; do
    if [ -f "$path" ]; then
        echo "✓ 找到适配器文件: $path"
        ADAPTER_FOUND=true
        break
    fi
done

if [ -z "$ADAPTER_FOUND" ]; then
    echo "⚠ 未找到Trimmomatic适配器文件，将手动下载..."
    
    # 创建适配器目录
    mkdir -p "$CONDA_PREFIX/share/trimmomatic/adapters"
    
    # 下载适配器文件
    wget https://github.com/timflutre/trimmomatic/raw/master/adapters/TruSeq3-PE-2.fa \
         -O "$CONDA_PREFIX/share/trimmomatic/adapters/TruSeq3-PE-2.fa"
    
    if [ -f "$CONDA_PREFIX/share/trimmomatic/adapters/TruSeq3-PE-2.fa" ]; then
        echo "✓ 适配器文件下载成功"
    else
        echo "❌ 适配器文件下载失败"
    fi
fi

echo ""
echo "=== 安装完成 ==="
echo "Conda环境: $ENV_NAME"
echo "激活环境: conda activate $ENV_NAME"
echo "退出环境: conda deactivate"
echo ""
echo "现在可以运行严格质控脚本了！"
