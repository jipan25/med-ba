#!/usr/bin/env python3
"""
质控优化执行器
实施质量修剪和自动化质控优化
"""

import subprocess
import os
from pathlib import Path
import pandas as pd
import sys

class QualityControlOptimizer:
    def __init__(self):
        self.data_dir = Path("../mnt-med/CRA007360")
        self.optimized_dir = Path("optimized_qc_results")
        self.optimized_dir.mkdir(exist_ok=True)
        
    def check_tools_availability(self):
        """检查必要工具是否可用"""
        print("=== 检查工具可用性 ===")
        
        required_tools = ['trim_galore', 'fastqc', 'multiqc']
        
        for tool in required_tools:
            try:
                result = subprocess.run([tool, '--version'], capture_output=True, text=True)
                if result.returncode == 0:
                    print(f"✓ {tool} 可用")
                else:
                    print(f"❌ {tool} 不可用，尝试安装...")
                    self.install_tool(tool)
            except FileNotFoundError:
                print(f"❌ {tool} 未安装，尝试安装...")
                self.install_tool(tool)
    
    def install_tool(self, tool):
        """安装必要的工具"""
        tool_install_commands = {
            'trim_galore': 'conda install -c bioconda trim-galore -y',
            'fastqc': 'conda install -c bioconda fastqc -y',
            'multiqc': 'pip install multiqc'
        }
        
        if tool in tool_install_commands:
            print(f"安装 {tool}...")
            try:
                subprocess.run(tool_install_commands[tool], shell=True, check=True)
                print(f"✓ {tool} 安装成功")
            except subprocess.CalledProcessError as e:
                print(f"❌ {tool} 安装失败: {e}")
                return False
        return True
    
    def discover_fastq_files(self):
        """发现FASTQ文件"""
        print("\n=== 发现FASTQ文件 ===")
        
        samples = {}
        
        # 查找样本目录
        sample_dirs = [d for d in self.data_dir.iterdir() if d.is_dir() and d.name.startswith("CRR")]
        
        for sample_dir in sample_dirs:
            sample_id = sample_dir.name
            
            # 查找FASTQ文件
            r1_files = list(sample_dir.glob("*_f1.fq.gz"))
            r2_files = list(sample_dir.glob("*_r2.fq.gz"))
            
            if r1_files and r2_files:
                samples[sample_id] = {
                    'r1': r1_files[0],
                    'r2': r2_files[0],
                    'group': 'BA' if '834' in sample_id else 'Control'
                }
                print(f"✓ {sample_id}: R1={r1_files[0].name}, R2={r2_files[0].name}")
        
        return samples
    
    def run_quality_trimming(self, samples):
        """运行质量修剪"""
        print("\n=== 执行质量修剪 ===")
        
        trimmed_samples = {}
        
        for sample_id, info in samples.items():
            print(f"处理样本 {sample_id}...")
            
            # 创建样本输出目录
            sample_out_dir = self.optimized_dir / sample_id
            sample_out_dir.mkdir(exist_ok=True)
            
            # trim_galore命令
            cmd = f"""
trim_galore \
    --quality 20 \
    --length 50 \
    --paired \
    --output_dir {sample_out_dir} \
    {info['r1']} \
    {info['r2']}
"""
            
            try:
                print(f"执行命令: trim_galore --quality 20 --length 50 --paired {info['r1'].name} {info['r2'].name}")
                result = subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True)
                
                # 检查修剪后的文件
                trimmed_r1 = sample_out_dir / f"{info['r1'].stem}_val_1.fq.gz"
                trimmed_r2 = sample_out_dir / f"{info['r2'].stem}_val_2.fq.gz"
                
                if trimmed_r1.exists() and trimmed_r2.exists():
                    trimmed_samples[sample_id] = {
                        'r1': trimmed_r1,
                        'r2': trimmed_r2,
                        'group': info['group'],
                        'trimming_report': sample_out_dir / f"{info['r1'].stem}_trimming_report.txt"
                    }
                    print(f"✓ {sample_id} 修剪完成")
                else:
                    print(f"❌ {sample_id} 修剪失败，输出文件不存在")
                    
            except subprocess.CalledProcessError as e:
                print(f"❌ {sample_id} 修剪失败: {e}")
                print(f"错误输出: {e.stderr}")
        
        return trimmed_samples
    
    def run_quality_control(self, trimmed_samples):
        """运行质控验证"""
        print("\n=== 运行质控验证 ===")
        
        qc_dir = self.optimized_dir / "fastqc_reports"
        qc_dir.mkdir(exist_ok=True)
        
        for sample_id, info in trimmed_samples.items():
            print(f"质控验证 {sample_id}...")
            
            cmd = f"fastqc {info['r1']} {info['r2']} -o {qc_dir}"
            
            try:
                result = subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True)
                print(f"✓ {sample_id} 质控完成")
            except subprocess.CalledProcessError as e:
                print(f"❌ {sample_id} 质控失败: {e}")
    
    def generate_multiqc_report(self):
        """生成MultiQC汇总报告"""
        print("\n=== 生成MultiQC汇总报告 ===")
        
        qc_dir = self.optimized_dir / "fastqc_reports"
        
        if not any(qc_dir.glob("*.html")):
            print("❌ 未找到质控报告，跳过MultiQC")
            return
        
        cmd = f"multiqc {qc_dir} -o {self.optimized_dir}"
        
        try:
            result = subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True)
            print("✓ MultiQC报告生成完成")
        except subprocess.CalledProcessError as e:
            print(f"❌ MultiQC失败: {e}")
    
    def validate_optimization_results(self):
        """验证优化结果"""
        print("\n=== 验证优化结果 ===")
        
        # 检查修剪后的文件
        trimmed_files = list(self.optimized_dir.rglob("*_val_*.fq.gz"))
        print(f"修剪后文件数量: {len(trimmed_files)}")
        
        # 检查质控报告
        qc_reports = list(self.optimized_dir.rglob("*fastqc.html"))
        print(f"质控报告数量: {len(qc_reports)}")
        
        # 检查MultiQC报告
        multiqc_report = self.optimized_dir / "multiqc_report.html"
        if multiqc_report.exists():
            print("✓ MultiQC汇总报告已生成")
        
        return len(trimmed_files) > 0 and len(qc_reports) > 0
    
    def generate_optimization_report(self, samples, trimmed_samples):
        """生成优化报告"""
        print("\n=== 生成优化报告 ===")
        
        report_content = f"""
# 质控优化执行报告

## 执行摘要
- 原始样本数: {len(samples)}
- 成功修剪样本数: {len(trimmed_samples)}
- 优化完成时间: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}

## 优化参数
- 质量阈值: 20 (Phred分数)
- 最小长度: 50 bp
- 接头去除: 自动检测
- 输出格式: 压缩FASTQ

## 样本处理状态

### 成功处理的样本
{chr(10).join(f"- {sample_id} ({info['group']})" for sample_id, info in trimmed_samples.items())}

### 处理失败的样本
{chr(10).join(f"- {sample_id}" for sample_id in set(samples.keys()) - set(trimmed_samples.keys())) if set(samples.keys()) - set(trimmed_samples.keys()) else "无"}

## 文件输出
- 修剪后文件: {self.optimized_dir}/
- 质控报告: {self.optimized_dir}/fastqc_reports/
- 汇总报告: {self.optimized_dir}/multiqc_report.html

## 下一步建议
1. 检查MultiQC报告确认质控通过
2. 进行序列比对分析
3. 执行基因表达定量
"""
        
        with open(self.optimized_dir / "optimization_report.md", "w") as f:
            f.write(report_content)
        
        print(f"✓ 优化报告已生成: {self.optimized_dir / 'optimization_report.md'}")

def main():
    """主执行流程"""
    print("=== 质控优化执行开始 ===")
    
    optimizer = QualityControlOptimizer()
    
    try:
        # 1. 检查工具
        optimizer.check_tools_availability()
        
        # 2. 发现文件
        samples = optimizer.discover_fastq_files()
        if not samples:
            print("❌ 未发现FASTQ文件")
            return
        
        # 3. 质量修剪
        trimmed_samples = optimizer.run_quality_trimming(samples)
        if not trimmed_samples:
            print("❌ 所有样本修剪失败")
            return
        
        # 4. 质控验证
        optimizer.run_quality_control(trimmed_samples)
        
        # 5. 生成汇总报告
        optimizer.generate_multiqc_report()
        
        # 6. 验证结果
        success = optimizer.validate_optimization_results()
        
        # 7. 生成报告
        optimizer.generate_optimization_report(samples, trimmed_samples)
        
        if success:
            print("\n🎉 质控优化执行成功!")
            print(f"所有结果保存在: {optimizer.optimized_dir}")
        else:
            print("\n⚠️ 质控优化部分完成，请检查结果")
            
    except Exception as e:
        print(f"❌ 质控优化执行失败: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()