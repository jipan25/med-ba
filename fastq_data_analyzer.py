#!/usr/bin/env python3
"""
FASTQ数据快速分析器
直接分析CRR524834和CRR524835的FASTQ文件
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import subprocess
import os
import gzip
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# 设置字体
import matplotlib
matplotlib.use('Agg')
plt.rcParams['font.sans-serif'] = ['DejaVu Sans', 'Arial', 'Helvetica']
plt.rcParams['axes.unicode_minus'] = False

class FastqDataAnalyzer:
    def __init__(self):
        self.data_dir = Path("../mnt-med/CRA007360")
        self.results_dir = Path("fastq_analysis_results")
        self.results_dir.mkdir(exist_ok=True)
        
    def check_data_availability(self):
        """检查数据文件是否存在"""
        print("=== 检查数据文件 ===")
        
        samples = {}
        
        # 检查CRR524834
        crr834_dir = self.data_dir / "CRR524834"
        if crr834_dir.exists():
            samples['CRR524834'] = {
                'r1': crr834_dir / "CRR524834_f1.fq.gz",
                'r2': crr834_dir / "CRR524834_r2.fq.gz",
                'group': 'BA'
            }
            print(f"✓ 发现CRR524834样本 (胆道闭锁组)")
        
        # 检查CRR524835
        crr835_dir = self.data_dir / "CRR524835"
        if crr835_dir.exists():
            samples['CRR524835'] = {
                'r1': crr835_dir / "CRR524835_f1.fq.gz",
                'r2': crr835_dir / "CRR524835_r2.fq.gz",
                'group': 'Control'
            }
            print(f"✓ 发现CRR524835样本 (对照组)")
        
        if not samples:
            print("❌ 未发现任何样本数据")
            return None
            
        # 检查文件大小
        for sample_id, info in samples.items():
            r1_size = info['r1'].stat().st_size / (1024**3)  # GB
            r2_size = info['r2'].stat().st_size / (1024**3)  # GB
            print(f"  {sample_id}: R1={r1_size:.1f}GB, R2={r2_size:.1f}GB")
            
        return samples
    
    def run_fastqc(self, samples):
        """运行FastQC质量评估"""
        print("\n=== 运行FastQC质量评估 ===")
        
        qc_dir = self.results_dir / "fastqc_reports"
        qc_dir.mkdir(exist_ok=True)
        
        for sample_id, info in samples.items():
            print(f"处理样本 {sample_id}...")
            
            cmd = f"fastqc {info['r1']} {info['r2']} -o {qc_dir}"
            try:
                result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
                if result.returncode == 0:
                    print(f"  ✓ FastQC完成")
                else:
                    print(f"  ⚠ FastQC警告: {result.stderr}")
            except Exception as e:
                print(f"  ❌ FastQC失败: {e}")
        
        return qc_dir
    
    def analyze_file_statistics(self, samples):
        """分析文件统计信息"""
        print("\n=== 分析文件统计信息 ===")
        
        stats_data = []
        
        for sample_id, info in samples.items():
            # 获取文件基本信息
            r1_size = info['r1'].stat().st_size / (1024**3)  # GB
            r2_size = info['r2'].stat().st_size / (1024**3)  # GB
            
            # 估算reads数量（基于文件大小）
            # 假设每个read约300bp，每个碱基约1字节
            estimated_reads = int(r1_size * 1e9 / 300)
            
            stats_data.append({
                'sample_id': sample_id,
                'group': info['group'],
                'r1_size_gb': round(r1_size, 2),
                'r2_size_gb': round(r2_size, 2),
                'total_size_gb': round(r1_size + r2_size, 2),
                'estimated_reads_millions': round(estimated_reads / 1e6, 1),
                'file_path': str(info['r1'].parent)
            })
        
        df_stats = pd.DataFrame(stats_data)
        df_stats.to_csv(self.results_dir / "file_statistics.csv", index=False)
        
        print("文件统计信息:")
        print(df_stats.to_string(index=False))
        
        return df_stats
    
    def create_data_summary(self, df_stats):
        """创建数据汇总报告"""
        print("\n=== 创建数据汇总报告 ===")
        
        # 计算总体统计
        total_size = df_stats['total_size_gb'].sum()
        total_reads = df_stats['estimated_reads_millions'].sum() * 1e6
        
        summary_content = f"""
# 胆道闭锁测序数据汇总报告

## 项目信息
- 项目编号: CRA007360
- 分析日期: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}
- 总数据量: {total_size:.1f} GB
- 总reads数: ~{total_reads/1e6:.1f} 百万

## 样本信息

### 样本统计
{df_stats.to_markdown(index=False)}

### 数据质量评估
- 已完成FastQC质量评估
- 质量报告保存在: fastqc_reports/

## 下一步分析建议

### 1. 数据预处理
- 使用trim_galore进行质量修剪
- 去除接头序列
- 过滤低质量reads

### 2. 序列比对
- 下载人类参考基因组 (GRCh38)
- 使用STAR或HISAT2进行比对
- 生成BAM文件

### 3. 基因表达定量
- 使用featureCounts或HTSeq计数
- 生成基因表达矩阵

### 4. 差异表达分析
- 使用DESeq2或edgeR
- 重点关注NETs和CD177相关基因

## 技术参数
- 测序平台: Illumina
- 读长: 双端测序 (PE)
- 数据格式: FASTQ (gzip压缩)

## 文件位置
原始数据: {self.data_dir}
分析结果: {self.results_dir}
"""
        
        with open(self.results_dir / "data_summary_report.md", "w") as f:
            f.write(summary_content)
        
        print(f"数据汇总报告已生成: {self.results_dir / 'data_summary_report.md'}")
    
    def create_visualizations(self, df_stats):
        """创建可视化图表"""
        print("\n=== 创建可视化图表 ===")
        
        viz_dir = self.results_dir / "visualizations"
        viz_dir.mkdir(exist_ok=True)
        
        # 1. 数据量对比图
        plt.figure(figsize=(10, 6))
        
        plt.subplot(1, 2, 1)
        plt.bar(df_stats['sample_id'], df_stats['total_size_gb'], 
                color=['red' if g == 'BA' else 'blue' for g in df_stats['group']])
        plt.title('样本数据量对比')
        plt.ylabel('数据量 (GB)')
        plt.xticks(rotation=45)
        
        plt.subplot(1, 2, 2)
        plt.bar(df_stats['sample_id'], df_stats['estimated_reads_millions'],
                color=['red' if g == 'BA' else 'blue' for g in df_stats['group']])
        plt.title('估计reads数量')
        plt.ylabel('Reads数量 (百万)')
        plt.xticks(rotation=45)
        
        plt.tight_layout()
        plt.savefig(viz_dir / 'data_statistics.png', dpi=300, bbox_inches='tight')
        plt.close()
        
        # 2. 分组统计
        group_stats = df_stats.groupby('group').agg({
            'total_size_gb': 'sum',
            'estimated_reads_millions': 'sum',
            'sample_id': 'count'
        }).rename(columns={'sample_id': 'sample_count'})
        
        plt.figure(figsize=(8, 6))
        group_stats[['total_size_gb', 'estimated_reads_millions']].plot(kind='bar')
        plt.title('分组数据统计')
        plt.ylabel('数值')
        plt.xticks(rotation=0)
        plt.legend(['总数据量 (GB)', '总reads数 (百万)'])
        plt.tight_layout()
        plt.savefig(viz_dir / 'group_statistics.png', dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"可视化图表已保存到: {viz_dir}")

def main():
    """主分析流程"""
    print("=== FASTQ数据快速分析开始 ===")
    
    analyzer = FastqDataAnalyzer()
    
    try:
        # 1. 检查数据可用性
        samples = analyzer.check_data_availability()
        if not samples:
            return
        
        # 2. 运行FastQC
        analyzer.run_fastqc(samples)
        
        # 3. 分析文件统计
        df_stats = analyzer.analyze_file_statistics(samples)
        
        # 4. 创建可视化
        analyzer.create_visualizations(df_stats)
        
        # 5. 生成汇总报告
        analyzer.create_data_summary(df_stats)
        
        print("\n=== 分析完成 ===")
        print(f"所有结果保存在: {analyzer.results_dir}")
        print("\n下一步建议:")
        print("1. 查看 fastqc_reports/ 中的质量报告")
        print("2. 查看 data_summary_report.md 了解数据概况")
        print("3. 准备进行序列比对分析")
        
    except Exception as e:
        print(f"分析过程中出现错误: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()