#!/bin/bash

# 严格质控优化管道（包含Conda环境设置）
# 针对低质量起始位置数据进行严格质控处理

set -euo pipefail

echo "=== 严格质控优化管道启动 ==="

# 配置参数
QUALITY_THRESHOLD=25
MIN_LENGTH=75
LEADING_BASES=15
TRAILING_BASES=10
SLIDING_WINDOW=4:20

DATA_DIR="../data/med-data/CRA007360"
OUTPUT_DIR="../data/results/strict_qc_results_$(date +%Y%m%d_%H%M%S)"

# Conda环境配置
CONDA_ENV_NAME="strict_qc_env"
CONDA_CHANNELS=("bioconda" "conda-forge")

echo "严格配置参数:"
echo "- 质量阈值: Q$QUALITY_THRESHOLD"
echo "- 最小长度: $MIN_LENGTH bp"
echo "- 起始修剪: $LEADING_BASES bp"
echo "- 滑动窗口: $SLIDING_WINDOW"
echo "- Conda环境: $CONDA_ENV_NAME"
echo "- 数据目录: $DATA_DIR"
echo "- 输出目录: $OUTPUT_DIR"

# 检查Conda是否安装
echo ""
echo "=== 检查Conda环境 ==="

check_conda() {
    if command -v conda &> /dev/null; then
        echo "✓ Conda 已安装"
        CONDA_AVAILABLE=true
        return 0
    else
        echo "❌ Conda 未安装"
        echo "请先安装Conda: https://docs.conda.io/en/latest/miniconda.html"
        exit 1
    fi
}

check_conda

# 配置Conda通道
setup_conda_channels() {
    echo "配置Conda通道..."
    for channel in "${CONDA_CHANNELS[@]}"; do
        if conda config --show channels | grep -q "$channel"; then
            echo "✓ 通道 $channel 已配置"
        else
            echo "添加通道: $channel"
            conda config --add channels "$channel"
        fi
    done
    conda config --set channel_priority strict
}

setup_conda_channels

# 创建或激活Conda环境
setup_conda_environment() {
    echo "设置Conda环境: $CONDA_ENV_NAME"
    
    if conda env list | grep -q "$CONDA_ENV_NAME"; then
        echo "✓ 环境 $CONDA_ENV_NAME 已存在，激活环境"
        # 注意：在脚本中激活conda环境需要使用conda activate
        eval "$(conda shell.bash hook)"
        conda activate "$CONDA_ENV_NAME"
    else
        echo "创建新环境: $CONDA_ENV_NAME"
        conda create -n "$CONDA_ENV_NAME" -y python=3.8
        eval "$(conda shell.bash hook)"
        conda activate "$CONDA_ENV_NAME"
        
        # 安装必要的工具
        echo "安装生物信息学工具..."
        conda install -y -c bioconda \
            trimmomatic=0.39 \
            fastqc=0.11.9 \
            multiqc=1.11 \
            fastp=0.23.2 \
            cutadapt=4.1
        
        pip install multiqc
        echo "✓ 所有工具安装完成"
    fi
}

setup_conda_environment

# 创建输出目录结构
mkdir -p "$OUTPUT_DIR"
mkdir -p "$OUTPUT_DIR/logs"
mkdir -p "$OUTPUT_DIR/qc_reports"

# 检查工具可用性
echo ""
echo "=== 检查工具可用性 ==="

check_tool() {
    if command -v "$1" &> /dev/null; then
        echo "✓ $1 可用"
        return 0
    else
        echo "❌ $1 未安装"
        return 1
    fi
}

for tool in trimmomatic fastqc multiqc fastp; do
    if ! check_tool "$tool"; then
        echo "正在安装 $tool..."
        conda install -c bioconda "$tool" -y
    fi
done

# 发现FASTQ文件
echo ""
echo "=== 发现FASTQ文件 ==="

declare -a SAMPLE_IDS=()
declare -A SAMPLE_FILES=()

for sample_dir in "$DATA_DIR"/CRR*; do
    if [ -d "$sample_dir" ]; then
        sample_id=$(basename "$sample_dir")
        r1_file="$sample_dir/${sample_id}_f1.fq.gz"
        r2_file="$sample_dir/${sample_id}_r2.fq.gz"

        if [ -f "$r1_file" ] && [ -f "$r2_file" ]; then
            SAMPLE_IDS+=("$sample_id")
            SAMPLE_FILES["${sample_id}_R1"]="$r1_file"
            SAMPLE_FILES["${sample_id}_R2"]="$r2_file"
            echo "✓ 发现样本: $sample_id"
            
            # 显示文件信息
            r1_size=$(du -h "$r1_file" | cut -f1)
            r2_size=$(du -h "$r2_file" | cut -f1)
            echo "  文件: $(basename "$r1_file") (${r1_size}), $(basename "$r2_file") (${r2_size})"
        fi
    fi
done

if [ ${#SAMPLE_IDS[@]} -eq 0 ]; then
    echo "❌ 未发现任何样本"
    exit 1
fi

echo "总共发现 ${#SAMPLE_IDS[@]} 个样本"

# 文件完整性检查函数
check_file_integrity() {
    local file="$1"
    local log_file="$2"
    
    if [ ! -f "$file" ]; then
        echo "文件不存在: $file" >> "$log_file"
        return 1
    fi

    local size
    size=$(stat -c%s "$file" 2>/dev/null || stat -f%z "$file" 2>/dev/null)
    
    if [ "$size" -lt 1024 ]; then
        echo "文件大小异常: $file (大小: ${size} bytes)" >> "$log_file"
        return 1
    fi

    if [[ "$file" == *.gz ]]; then
        if ! gzip -t "$file" 2>> "$log_file"; then
            echo "gzip文件损坏: $file" >> "$log_file"
            return 1
        fi
    fi

    echo "文件检查通过: $file (大小: ${size} bytes)" >> "$log_file"
    return 0
}

# 第一步：原始数据质量评估
echo ""
echo "=== 第一步：原始数据质量评估 ==="

RAW_QC_DIR="$OUTPUT_DIR/raw_fastqc"
mkdir -p "$RAW_QC_DIR"

for sample_id in "${SAMPLE_IDS[@]}"; do
    echo "评估原始数据质量: $sample_id"
    
    r1_file="${SAMPLE_FILES["${sample_id}_R1"]}"
    r2_file="${SAMPLE_FILES["${sample_id}_R2"]}"
    
    log_file="$OUTPUT_DIR/logs/${sample_id}_raw_qc.log"
    
    # 检查文件完整性
    if ! check_file_integrity "$r1_file" "$log_file" || ! check_file_integrity "$r2_file" "$log_file"; then
        echo "❌ $sample_id 原始文件完整性检查失败"
        continue
    fi
    
    # 运行FastQC
    if fastqc "$r1_file" "$r2_file" -o "$RAW_QC_DIR" -t 2 --quiet 2>> "$log_file"; then
        echo "✓ $sample_id 原始数据质量评估完成"
    else
        echo "❌ $sample_id 原始数据质量评估失败"
    fi
done

# 第二步：多重质量修剪策略
echo ""
echo "=== 第二步：多重质量修剪策略 ==="

# 查找Trimmomatic适配器文件
find_adapters() {
    local possible_paths=(
        "$(dirname $(which trimmomatic))/../share/trimmomatic/adapters/TruSeq3-PE-2.fa"
        "$CONDA_PREFIX/share/trimmomatic/adapters/TruSeq3-PE-2.fa"
        "/usr/share/trimmomatic/adapters/TruSeq3-PE-2.fa"
        "./adapters/TruSeq3-PE-2.fa"
    )
    
    for path in "${possible_paths[@]}"; do
        if [ -f "$path" ]; then
            echo "$path"
            return 0
        fi
    done
    
    echo "未找到适配器文件，将使用HEADCROP代替ILLUMINACLIP" >&2
    echo ""
    return 1
}

ADAPTER_FILE=$(find_adapters)

for sample_id in "${SAMPLE_IDS[@]}"; do
    echo "处理样本: $sample_id"
    
    r1_file="${SAMPLE_FILES["${sample_id}_R1"]}"
    r2_file="${SAMPLE_FILES["${sample_id}_R2"]}"
    
    sample_out_dir="$OUTPUT_DIR/$sample_id"
    mkdir -p "$sample_out_dir"
    
    log_file="$OUTPUT_DIR/logs/${sample_id}_trimming.log"
    
    # 输出文件定义
    trimmed_r1_paired="$sample_out_dir/${sample_id}_R1_trimmed_paired.fq.gz"
    trimmed_r1_unpaired="$sample_out_dir/${sample_id}_R1_trimmed_unpaired.fq.gz"
    trimmed_r2_paired="$sample_out_dir/${sample_id}_R2_trimmed_paired.fq.gz"
    trimmed_r2_unpaired="$sample_out_dir/${sample_id}_R2_trimmed_unpaired.fq.gz"
    
    # 检查是否已完成修剪
    if [ -f "$trimmed_r1_paired" ] && [ -f "$trimmed_r2_paired" ]; then
        if check_file_integrity "$trimmed_r1_paired" "$log_file" && check_file_integrity "$trimmed_r2_paired" "$log_file"; then
            echo "✓ $sample_id 已经完成修剪，跳过"
            continue
        else
            echo "⚠ $sample_id 修剪文件不完整，重新处理"
            rm -f "$trimmed_r1_paired" "$trimmed_r2_paired" "$trimmed_r1_unpaired" "$trimmed_r2_unpaired"
        fi
    fi
    
    # 构建Trimmomatic命令
    if [ -n "$ADAPTER_FILE" ] && [ -f "$ADAPTER_FILE" ]; then
        echo "使用适配器文件: $ADAPTER_FILE"
        trimmomatic_cmd="trimmomatic PE -threads 12 -phred33 \
            '$r1_file' '$r2_file' \
            '$trimmed_r1_paired' '$trimmed_r1_unpaired' \
            '$trimmed_r2_paired' '$trimmed_r2_unpaired' \
            ILLUMINACLIP:${ADAPTER_FILE}:2:30:10:8:true \
            HEADCROP:$LEADING_BASES \
            LEADING:$QUALITY_THRESHOLD \
            TRAILING:$QUALITY_THRESHOLD \
            SLIDINGWINDOW:$SLIDING_WINDOW \
            MINLEN:$MIN_LENGTH"
    else
        echo "未找到适配器文件，跳过接头修剪"
        trimmomatic_cmd="trimmomatic PE -threads 12 -phred33 \
            '$r1_file' '$r2_file' \
            '$trimmed_r1_paired' '$trimmed_r1_unpaired' \
            '$trimmed_r2_paired' '$trimmed_r2_unpaired' \
            HEADCROP:$LEADING_BASES \
            LEADING:$QUALITY_THRESHOLD \
            TRAILING:$QUALITY_THRESHOLD \
            SLIDINGWINDOW:$SLIDING_WINDOW \
            MINLEN:$MIN_LENGTH"
    fi
    
    echo "执行严格质量修剪: $sample_id"
    echo "执行命令: $trimmomatic_cmd" >> "$log_file"
    
    if eval "$trimmomatic_cmd" >> "$log_file" 2>&1; then
        # 检查修剪后文件
        if check_file_integrity "$trimmed_r1_paired" "$log_file" && check_file_integrity "$trimmed_r2_paired" "$log_file"; then
            # 统计修剪结果
            original_reads=$(zcat "$r1_file" 2>/dev/null | awk 'END {print NR/4}')
            trimmed_reads=$(zcat "$trimmed_r1_paired" 2>/dev/null | awk 'END {print NR/4}')
            
            if [ "$original_reads" -gt 0 ]; then
                survival_rate=$(echo "scale=2; $trimmed_reads * 100 / $original_reads" | bc)
                echo "✓ $sample_id 修剪完成 (存活率: ${survival_rate}%)" | tee -a "$log_file"
            else
                echo "✓ $sample_id 修剪完成" | tee -a "$log_file"
            fi
            echo "原始reads: $original_reads, 修剪后reads: $trimmed_reads" >> "$log_file"
        else
            echo "❌ $sample_id 修剪后文件检查失败" | tee -a "$log_file"
        fi
    else
        echo "❌ $sample_id 修剪失败" | tee -a "$log_file"
    fi
done

# 第三步：修剪后质量验证
echo ""
echo "=== 第三步：修剪后质量验证 ==="

TRIMMED_QC_DIR="$OUTPUT_DIR/trimmed_fastqc"
mkdir -p "$TRIMMED_QC_DIR"

declare -a SUCCESSFUL_SAMPLES=()

for sample_id in "${SAMPLE_IDS[@]}"; do
    sample_out_dir="$OUTPUT_DIR/$sample_id"
    trimmed_r1_paired="$sample_out_dir/${sample_id}_R1_trimmed_paired.fq.gz"
    trimmed_r2_paired="$sample_out_dir/${sample_id}_R2_trimmed_paired.fq.gz"
    
    log_file="$OUTPUT_DIR/logs/${sample_id}_trimmed_qc.log"
    
    # 检查修剪文件是否存在且完整
    if [ ! -f "$trimmed_r1_paired" ] || [ ! -f "$trimmed_r2_paired" ]; then
        echo "❌ $sample_id 修剪文件不存在，跳过质控"
        continue
    fi
    
    if ! check_file_integrity "$trimmed_r1_paired" "$log_file" || ! check_file_integrity "$trimmed_r2_paired" "$log_file"; then
        echo "❌ $sample_id 修剪文件不完整，跳过质控"
        continue
    fi
    
    echo "质控验证: $sample_id"
    
    if fastqc "$trimmed_r1_paired" "$trimmed_r2_paired" -o "$TRIMMED_QC_DIR" -t 2 --quiet 2>> "$log_file"; then
        echo "✓ $sample_id 修剪后质控完成"
        SUCCESSFUL_SAMPLES+=("$sample_id")
    else
        echo "❌ $sample_id 修剪后质控失败"
    fi
done

# 第四步：生成质量报告
echo ""
echo "=== 第四步：生成质量报告 ==="

# 生成原始数据MultiQC报告
if [ -d "$RAW_QC_DIR" ] && [ "$(ls -A "$RAW_QC_DIR")" ]; then
    echo "生成原始数据MultiQC报告..."
    multiqc "$RAW_QC_DIR" -o "$OUTPUT_DIR" -n "raw_multiqc_report.html" --quiet
fi

# 生成修剪后数据MultiQC报告
if [ -d "$TRIMMED_QC_DIR" ] && [ "$(ls -A "$TRIMMED_QC_DIR")" ]; then
    echo "生成修剪后数据MultiQC报告..."
    multiqc "$TRIMMED_QC_DIR" -o "$OUTPUT_DIR" -n "trimmed_multiqc_report.html" --quiet
fi

# 第五步：生成详细执行报告
echo ""
echo "=== 第五步：生成详细执行报告 ==="

REPORT_FILE="$OUTPUT_DIR/strict_qc_execution_report.md"

{
echo "# 严格质控优化管道执行报告"
echo ""
echo "## 执行摘要"
echo "- 执行时间: $(date)"
echo "- 总样本数: ${#SAMPLE_IDS[@]}"
echo "- 成功处理样本数: ${#SUCCESSFUL_SAMPLES[@]}"
success_rate=$(( ${#SUCCESSFUL_SAMPLES[@]} * 100 / ${#SAMPLE_IDS[@]} ))
echo "- 成功率: ${success_rate}%"
echo "- Conda环境: $CONDA_ENV_NAME"
echo ""
echo "## 严格质控参数"
echo "- 质量阈值: Q$QUALITY_THRESHOLD"
echo "- 最小长度: $MIN_LENGTH bp"
echo "- 起始修剪: $LEADING_BASES bp"
echo "- 滑动窗口: $SLIDING_WINDOW"
echo "- 适配器文件: ${ADAPTER_FILE:-未使用}"
echo ""
echo "## 处理结果统计"
echo ""
echo "### 成功处理的样本"
for sample_id in "${SUCCESSFUL_SAMPLES[@]}"; do
    sample_out_dir="$OUTPUT_DIR/$sample_id"
    trimmed_r1_paired="$sample_out_dir/${sample_id}_R1_trimmed_paired.fq.gz"
    
    if [ -f "$trimmed_r1_paired" ]; then
        read_count=$(zcat "$trimmed_r1_paired" 2>/dev/null | awk 'END {print NR/4}')
        echo "- $sample_id (reads: $read_count)"
    else
        echo "- $sample_id"
    fi
done
echo ""
echo "### 处理失败的样本"
for sample_id in "${SAMPLE_IDS[@]}"; do
    if [[ ! " ${SUCCESSFUL_SAMPLES[@]} " =~ " ${sample_id} " ]]; then
        echo "- $sample_id"
    fi
done
echo ""
echo "## 环境信息"
echo "\`\`\`"
conda info >> "$REPORT_FILE.tmp" 2>&1
cat "$REPORT_FILE.tmp"
rm "$REPORT_FILE.tmp"
echo "\`\`\`"
echo ""
echo "## 安装的工具版本"
echo "\`\`\`"
for tool in trimmomatic fastqc multiqc fastp; do
    if command -v "$tool" &> /dev/null; then
        echo -n "$tool: "
        $tool --version 2>/dev/null | head -1 || echo "版本信息不可用"
    fi
done
echo "\`\`\`"
echo ""
echo "## 下一步建议"
echo "1. 检查 \`trimmed_multiqc_report.html\` 确认质控效果"
echo "2. 特别关注 per-base quality 图表，确认起始位置质量已改善"
echo "3. 如需重新运行，可删除输出目录或调整参数："
echo "   - 提高质量阈值至Q30"
echo "   - 增加起始修剪长度至20"
echo "   - 使用更严格的滑动窗口(4:25)"
echo "4. 质控通过后进行序列比对分析"
} > "$REPORT_FILE"

echo "✓ 详细执行报告已生成: $REPORT_FILE"

# 最终统计和验证
echo ""
echo "=== 最终验证 ==="

echo "原始数据质控报告: $(find "$RAW_QC_DIR" -name "*.html" 2>/dev/null | wc -l)"
echo "修剪后质控报告: $(find "$TRIMMED_QC_DIR" -name "*.html" 2>/dev/null | wc -l)"
echo "成功修剪样本: ${#SUCCESSFUL_SAMPLES[@]}/${#SAMPLE_IDS[@]}"

if [ -f "$OUTPUT_DIR/trimmed_multiqc_report.html" ]; then
    echo "✓ 修剪后MultiQC报告存在"
else
    echo "⚠ 修剪后MultiQC报告未生成"
fi

echo ""
echo "🎉 严格质控优化管道执行完成!"
echo ""
echo "📊 重要提醒:"
echo "   请打开 $OUTPUT_DIR/trimmed_multiqc_report.html"
echo "   重点检查 'Per base sequence quality' 图表"
echo "   确认起始位置质量已得到显著改善"
echo ""
echo "🔧 环境管理:"
echo "   Conda环境: $CONDA_ENV_NAME"
echo "   激活环境: conda activate $CONDA_ENV_NAME"
echo "   退出环境: conda deactivate"
echo ""
echo "输出目录: $OUTPUT_DIR"

# 停用conda环境
conda deactivate
