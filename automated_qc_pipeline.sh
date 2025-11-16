#!/bin/bash

# 自动化质控优化管道
# 执行完整的质控优化流程，支持断点续执行

set -e  # 遇到错误立即退出

# 确保使用bash
if [ -z "$BASH_VERSION" ]; then
    echo "错误: 请使用bash运行此脚本"
    exit 1
fi

echo "=== 自动化质控优化管道启动 ==="

# 配置参数
QUALITY_THRESHOLD=20
MIN_LENGTH=50
DATA_DIR="../mnt-med/CRA007360"
OUTPUT_DIR="../mnt-med-data/optimized_qc_results_20251115_195620"  # 使用固定目录名以便续执行

echo "配置参数:"
echo "- 质量阈值: $QUALITY_THRESHOLD"
echo "- 最小长度: $MIN_LENGTH"
echo "- 数据目录: $DATA_DIR"
echo "- 输出目录: $OUTPUT_DIR"

# 创建输出目录
mkdir -p "$OUTPUT_DIR"

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

check_tool trim_galore || {
    echo "安装 trim_galore..."
    conda install -c bioconda trim-galore -y
}

check_tool fastqc || {
    echo "安装 fastqc..."
    conda install -c bioconda fastqc -y
}

check_tool multiqc || {
    echo "安装 multiqc..."
    pip install multiqc
}

# 发现FASTQ文件
echo ""
echo "=== 发现FASTQ文件 ==="

SAMPLES=""
SAMPLE_COUNT=0

for sample_dir in "$DATA_DIR"/CRR*; do
    if [ -d "$sample_dir" ]; then
        sample_id=$(basename "$sample_dir")
        r1_file="$sample_dir/${sample_id}_f1.fq.gz"
        r2_file="$sample_dir/${sample_id}_r2.fq.gz"
        
        if [ -f "$r1_file" ] && [ -f "$r2_file" ]; then
            if [ -z "$SAMPLES" ]; then
                SAMPLES="$sample_id"
            else
                SAMPLES="$SAMPLES $sample_id"
            fi
            SAMPLE_COUNT=$((SAMPLE_COUNT + 1))
            echo "✓ 发现样本: $sample_id"
        fi
    fi
done

if [ $SAMPLE_COUNT -eq 0 ]; then
    echo "❌ 未发现任何样本"
    exit 1
fi

echo "总共发现 $SAMPLE_COUNT 个样本"

# 检查文件完整性的函数
check_file_complete() {
    local file="$1"
    if [ ! -f "$file" ]; then
        return 1
    fi
    
    # 检查文件大小是否合理（至少1KB）
    local size=$(stat -c%s "$file" 2>/dev/null || stat -f%z "$file" 2>/dev/null)
    if [ "$size" -lt 1024 ]; then
        return 1
    fi
    
    # 检查gzip文件完整性
    if [[ "$file" == *.gz ]]; then
        if gzip -t "$file" 2>/dev/null; then
            return 0
        else
            return 1
        fi
    fi
    
    return 0
}

# 质量修剪（支持断点续执行）
echo ""
echo "=== 执行质量修剪 ==="

for sample_id in $SAMPLES; do
    echo "检查样本: $sample_id"
    
    r1_file="$DATA_DIR/$sample_id/${sample_id}_f1.fq.gz"
    r2_file="$DATA_DIR/$sample_id/${sample_id}_r2.fq.gz"
    
    sample_out_dir="$OUTPUT_DIR/$sample_id"
    mkdir -p "$sample_out_dir"
    
    trimmed_r1="$sample_out_dir/${sample_id}_f1_val_1.fq.gz"
    trimmed_r2="$sample_out_dir/${sample_id}_r2_val_2.fq.gz"
    
    # 检查是否已经完成修剪
    if check_file_complete "$trimmed_r1" && check_file_complete "$trimmed_r2"; then
        echo "✓ $sample_id 已经完成修剪，跳过"
        continue
    fi
    
    # 如果val文件存在但不完整，删除重新处理
    if [ -f "$trimmed_r1" ] && ! check_file_complete "$trimmed_r1"; then
        echo "⚠ $sample_id R1 val文件不完整，重新处理"
        rm -f "$trimmed_r1"
    fi
    
    if [ -f "$trimmed_r2" ] && ! check_file_complete "$trimmed_r2"; then
        echo "⚠ $sample_id R2 val文件不完整，重新处理"
        rm -f "$trimmed_r2"
    fi
    
    # 执行trim_galore
    echo "执行质量修剪: $sample_id"
    
    if trim_galore \
        --quality "$QUALITY_THRESHOLD" \
        --length "$MIN_LENGTH" \
        --paired \
        --cores 4 \
        --output_dir "$sample_out_dir" \
        "$r1_file" \
        "$r2_file"; then
        echo "✓ $sample_id 修剪完成"
    else
        echo "❌ $sample_id 修剪失败"
    fi
done

# 检查修剪结果
echo ""
echo "=== 检查修剪结果 ==="

TRIMMED_SAMPLES=""
TRIMMED_COUNT=0
for sample_id in $SAMPLES; do
    sample_out_dir="$OUTPUT_DIR/$sample_id"
    trimmed_r1="$sample_out_dir/${sample_id}_f1_val_1.fq.gz"
    trimmed_r2="$sample_out_dir/${sample_id}_r2_val_2.fq.gz"
    
    if check_file_complete "$trimmed_r1" && check_file_complete "$trimmed_r2"; then
        if [ -z "$TRIMMED_SAMPLES" ]; then
            TRIMMED_SAMPLES="$sample_id"
        else
            TRIMMED_SAMPLES="$TRIMMED_SAMPLES $sample_id"
        fi
        TRIMMED_COUNT=$((TRIMMED_COUNT + 1))
        echo "✓ $sample_id 修剪成功"
    else
        echo "❌ $sample_id 修剪失败或文件不完整"
    fi
done

if [ $TRIMMED_COUNT -eq 0 ]; then
    echo "❌ 所有样本修剪失败"
    exit 1
fi

# 质控验证（支持断点续执行）
echo ""
echo "=== 运行质控验证 ==="

QC_DIR="$OUTPUT_DIR/fastqc_reports"
mkdir -p "$QC_DIR"

for sample_id in $TRIMMED_SAMPLES; do
    sample_out_dir="$OUTPUT_DIR/$sample_id"
    trimmed_r1="$sample_out_dir/${sample_id}_f1_val_1.fq.gz"
    trimmed_r2="$sample_out_dir/${sample_id}_r2_val_2.fq.gz"
    
    # 检查是否已经生成了FastQC报告
    r1_report="$QC_DIR/$(basename "$trimmed_r1" .fq.gz)_fastqc.html"
    r2_report="$QC_DIR/$(basename "$trimmed_r2" .fq.gz)_fastqc.html"
    
    if [ -f "$r1_report" ] && [ -f "$r2_report" ]; then
        echo "✓ $sample_id 质控报告已存在，跳过"
        continue
    fi
    
    echo "质控验证: $sample_id"
    
    if fastqc "$trimmed_r1" "$trimmed_r2" -o "$QC_DIR" -t 4; then
        echo "✓ $sample_id 质控完成"
    else
        echo "❌ $sample_id 质控失败"
    fi
done

# 生成MultiQC报告（支持断点续执行）
echo ""
echo "=== 生成MultiQC汇总报告 ==="

MULTIQC_REPORT="$OUTPUT_DIR/multiqc_report.html"

if [ -f "$MULTIQC_REPORT" ]; then
    echo "✓ MultiQC报告已存在，跳过生成"
else
    if multiqc "$QC_DIR" -o "$OUTPUT_DIR"; then
        echo "✓ MultiQC报告生成完成"
    else
        echo "❌ MultiQC报告生成失败"
    fi
fi

# 生成执行报告
echo ""
echo "=== 生成执行报告 ==="

REPORT_FILE="$OUTPUT_DIR/pipeline_execution_report.md"

cat > "$REPORT_FILE" << EOF
# 自动化质控优化管道执行报告

## 执行摘要
- 执行时间: $(date)
- 原始样本数: $SAMPLE_COUNT
- 成功修剪样本数: $TRIMMED_COUNT
- 成功率: $((TRIMMED_COUNT * 100 / SAMPLE_COUNT))%

## 配置参数
- 质量阈值: $QUALITY_THRESHOLD
- 最小长度: $MIN_LENGTH
- 数据目录: $DATA_DIR
- 输出目录: $OUTPUT_DIR

## 样本处理状态

### 成功处理的样本
EOF

for sample_id in $TRIMMED_SAMPLES; do
    echo "- $sample_id" >> "$REPORT_FILE"
done

cat >> "$REPORT_FILE" << EOF

### 处理失败的样本
EOF

# 检查失败的样本
for sample_id in $SAMPLES; do
    found=0
    for trimmed_id in $TRIMMED_SAMPLES; do
        if [ "$sample_id" = "$trimmed_id" ]; then
            found=1
            break
        fi
    done
    if [ $found -eq 0 ]; then
        echo "- $sample_id" >> "$REPORT_FILE"
    fi
done

cat >> "$REPORT_FILE" << EOF

## 文件输出结构
- 修剪后FASTQ文件: $OUTPUT_DIR/{sample_id}/
- FastQC报告: $OUTPUT_DIR/fastqc_reports/
- MultiQC汇总报告: $OUTPUT_DIR/multiqc_report.html

## 下一步分析建议
1. 检查MultiQC报告确认质控通过
2. 进行序列比对分析
3. 执行基因表达定量
4. 差异表达分析

## 技术说明
- 使用trim_galore进行质量修剪和接头去除
- 使用FastQC进行质量评估
- 使用MultiQC生成汇总报告
- 管道设计为容错处理，单个样本失败不影响其他样本
EOF

echo "✓ 执行报告已生成: $REPORT_FILE"

# 最终验证
echo ""
echo "=== 最终验证 ==="

echo "修剪后文件数量: $(find "$OUTPUT_DIR" -name "*_val_*.fq.gz" | wc -l)"
echo "质控报告数量: $(find "$QC_DIR" -name "*.html" | wc -l)"

if [ -f "$OUTPUT_DIR/multiqc_report.html" ]; then
    echo "✓ MultiQC汇总报告存在"
else
    echo "❌ MultiQC汇总报告不存在"
fi

echo ""
echo "🎉 自动化质控优化管道执行完成!"
echo "输出目录: $OUTPUT_DIR"
echo ""
echo "下一步:"
echo "1. 打开 $OUTPUT_DIR/multiqc_report.html 查看质控结果"
echo "2. 确认质控通过后进行下游分析"
echo "3. 如有问题，查看 $REPORT_FILE"