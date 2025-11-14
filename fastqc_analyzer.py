#!/usr/bin/env python3
"""
FastQC报告解析器
专门解析质控失败原因并提供优化方案
"""

import re
import os
from pathlib import Path
from collections import defaultdict

class FastQCAnalyzer:
    def __init__(self, qc_dir):
        self.qc_dir = Path(qc_dir)
        self.results = {}
        
    def parse_fastqc_reports(self):
        """解析所有FastQC报告"""
        print("=== 解析FastQC质控报告 ===")
        
        html_files = list(self.qc_dir.glob("*.html"))
        
        for html_file in html_files:
            print(f"\n分析文件: {html_file.name}")
            sample_result = self.parse_single_report(html_file)
            self.results[html_file.name] = sample_result
            
        return self.results
    
    def parse_single_report(self, html_file):
        """解析单个FastQC报告"""
        with open(html_file, 'r', encoding='utf-8', errors='ignore') as f:
            content = f.read()
        
        # 提取质控模块状态
        module_status = self.extract_module_status(content)
        
        # 识别常见问题模式
        issues = self.identify_common_issues(content)
        
        return {
            'module_status': module_status,
            'issues': issues,
            'summary': self.generate_summary(module_status, issues)
        }
    
    def extract_module_status(self, content):
        """提取质控模块状态"""
        # 查找模块状态模式
        status_pattern = r'<td class="(FAIL|WARN|PASS)">([^<]+)</td>'
        matches = re.findall(status_pattern, content)
        
        module_status = {}
        for status, module_name in matches:
            # 清理模块名称
            module_name = module_name.strip()
            if module_name and len(module_name) > 2:  # 过滤掉空的和太短的
                module_status[module_name] = status
        
        return module_status
    
    def identify_common_issues(self, content):
        """识别常见质控问题"""
        issues = []
        
        # 检查常见问题关键词
        issue_patterns = {
            'low_quality': ['poor quality', 'quality drops', 'low scores'],
            'adapter_contamination': ['adapter content', 'illumina universal adapter'],
            'overrepresented_sequences': ['overrepresented sequences', 'overrepresented'],
            'kmer_content': ['kmer content', 'k-mer'],
            'per_tile_quality': ['per tile sequence quality', 'tile quality'],
            'sequence_length': ['sequence length distribution', 'length distribution']
        }
        
        for issue_type, keywords in issue_patterns.items():
            for keyword in keywords:
                if keyword.lower() in content.lower():
                    issues.append(issue_type)
                    break
        
        return list(set(issues))
    
    def generate_summary(self, module_status, issues):
        """生成质控摘要"""
        fail_modules = [mod for mod, status in module_status.items() if status == 'FAIL']
        warn_modules = [mod for mod, status in module_status.items() if status == 'WARN']
        pass_modules = [mod for mod, status in module_status.items() if status == 'PASS']
        
        summary = {
            'total_modules': len(module_status),
            'fail_count': len(fail_modules),
            'warn_count': len(warn_modules),
            'pass_count': len(pass_modules),
            'fail_modules': fail_modules,
            'warn_modules': warn_modules,
            'pass_modules': pass_modules,
            'issues_detected': issues
        }
        
        return summary
    
    def generate_optimization_report(self):
        """生成质控优化报告"""
        print("\n=== 质控优化分析报告 ===")
        
        # 汇总所有样本的问题
        all_fail_modules = set()
        all_issues = set()
        
        for filename, result in self.results.items():
            summary = result['summary']
            all_fail_modules.update(summary['fail_modules'])
            all_issues.update(summary['issues_detected'])
        
        # 生成优化建议
        optimization_suggestions = self.generate_suggestions(all_fail_modules, all_issues)
        
        report = f"""
# FastQC质控优化分析报告

## 质控结果汇总
- 分析样本数: {len(self.results)}
- 发现FAIL模块: {len(all_fail_modules)} 个
- 检测到问题类型: {len(all_issues)} 类

## 具体FAIL模块分析
{chr(10).join(f"- {module}" for module in all_fail_modules)}

## 检测到的问题类型
{chr(10).join(f"- {issue}" for issue in all_issues)}

## 优化建议

### 1. 数据修剪策略
{optimization_suggestions.get('trimming_strategy', '')}

### 2. 质控参数调整
{optimization_suggestions.get('qc_parameters', '')}

### 3. 下游分析注意事项
{optimization_suggestions.get('downstream_analysis', '')}

## 详细样本分析
"""
        
        # 添加每个样本的详细分析
        for filename, result in self.results.items():
            summary = result['summary']
            report += f"""
### {filename}
- FAIL模块: {len(summary['fail_modules'])} 个
- WARN模块: {len(summary['warn_modules'])} 个  
- PASS模块: {len(summary['pass_modules'])} 个
- 检测问题: {', '.join(summary['issues_detected']) if summary['issues_detected'] else '无'}
"""
        
        return report
    
    def generate_suggestions(self, fail_modules, issues):
        """根据问题生成优化建议"""
        suggestions = {
            'trimming_strategy': '',
            'qc_parameters': '',
            'downstream_analysis': ''
        }
        
        # 修剪策略建议
        trimming_suggestions = []
        if 'low_quality' in issues:
            trimming_suggestions.append("使用trim_galore进行质量修剪，设置--quality 20")
        if 'adapter_contamination' in issues:
            trimming_suggestions.append("启用自动接头检测和去除")
        if any('sequence_length' in mod.lower() for mod in fail_modules):
            trimming_suggestions.append("考虑设置最小读长过滤")
        
        suggestions['trimming_strategy'] = '\n'.join(trimming_suggestions) if trimming_suggestions else "当前数据质量良好，无需特殊修剪"
        
        # 质控参数建议
        qc_suggestions = []
        if 'overrepresented_sequences' in issues:
            qc_suggestions.append("检查是否存在污染序列，考虑使用kraken2进行物种鉴定")
        if 'kmer_content' in issues:
            qc_suggestions.append("注意k-mer偏差可能影响某些分析工具")
        
        suggestions['qc_parameters'] = '\n'.join(qc_suggestions) if qc_suggestions else "质控参数设置合理"
        
        # 下游分析建议
        downstream_suggestions = []
        if 'low_quality' in issues:
            downstream_suggestions.append("比对前建议进行质量修剪")
        if 'adapter_contamination' in issues:
            downstream_suggestions.append("确保比对前已去除接头序列")
        
        suggestions['downstream_analysis'] = '\n'.join(downstream_suggestions) if downstream_suggestions else "数据质量适合下游分析"
        
        return suggestions

def main():
    """主分析流程"""
    qc_dir = "analysis_results_1114/quality_control"
    
    if not os.path.exists(qc_dir):
        print(f"错误: 质控目录不存在: {qc_dir}")
        return
    
    analyzer = FastQCAnalyzer(qc_dir)
    results = analyzer.parse_fastqc_reports()
    
    # 生成优化报告
    report = analyzer.generate_optimization_report()
    
    # 保存报告
    output_file = "fastqc_optimization_report.md"
    with open(output_file, 'w', encoding='utf-8') as f:
        f.write(report)
    
    print(f"\n质控优化报告已生成: {output_file}")
    print("\n下一步:")
    print("1. 查看报告了解具体质控问题")
    print("2. 根据建议优化质控流程")
    print("3. 重新运行质控分析")

if __name__ == "__main__":
    main()