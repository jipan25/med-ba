#!/usr/bin/env python3
"""
详细FastQC报告解析器
直接解析HTML报告获取准确的质控状态
"""

import re
from pathlib import Path
from bs4 import BeautifulSoup
import os

def parse_fastqc_detailed(html_file):
    """详细解析单个FastQC报告"""
    print(f"\n=== 详细解析: {html_file.name} ===")
    
    with open(html_file, 'r', encoding='utf-8', errors='ignore') as f:
        content = f.read()
    
    # 使用BeautifulSoup解析HTML
    soup = BeautifulSoup(content, 'html.parser')
    
    # 查找所有模块状态
    modules = {}
    
    # 查找模块表格
    tables = soup.find_all('table')
    for table in tables:
        rows = table.find_all('tr')
        for row in rows:
            cells = row.find_all('td')
            if len(cells) >= 2:
                # 检查状态单元格
                status_cell = cells[0]
                module_cell = cells[1]
                
                # 获取状态类
                status_class = status_cell.get('class', [])
                if status_class:
                    status = status_class[0].upper()
                    module_name = module_cell.get_text(strip=True)
                    
                    if module_name and len(module_name) > 5:  # 过滤有效模块名
                        modules[module_name] = status
                        print(f"  {status}: {module_name}")
    
    return modules

def analyze_qc_issues(modules):
    """分析质控问题并提供优化建议"""
    fail_modules = [mod for mod, status in modules.items() if status == 'FAIL']
    warn_modules = [mod for mod, status in modules.items() if status == 'WARN']
    
    print(f"\n质控状态统计:")
    print(f"FAIL模块: {len(fail_modules)}")
    print(f"WARN模块: {len(warn_modules)}")
    print(f"PASS模块: {len(modules) - len(fail_modules) - len(warn_modules)}")
    
    # 分析具体问题
    issues = {}
    
    for module, status in modules.items():
        if status in ['FAIL', 'WARN']:
            module_lower = module.lower()
            
            if 'per base sequence quality' in module_lower:
                issues['base_quality'] = {
                    'severity': status,
                    'description': '碱基质量分数问题',
                    'solution': '使用trim_galore进行质量修剪，设置--quality 20'
                }
            
            elif 'per sequence quality scores' in module_lower:
                issues['sequence_quality'] = {
                    'severity': status,
                    'description': '序列质量分布问题',
                    'solution': '检查测序平台和文库构建质量'
                }
            
            elif 'per base sequence content' in module_lower:
                issues['base_content'] = {
                    'severity': status,
                    'description': '碱基组成偏差',
                    'solution': '可能是测序起始偏差，通常可接受'
                }
            
            elif 'adapter content' in module_lower:
                issues['adapter'] = {
                    'severity': status,
                    'description': '接头污染',
                    'solution': '使用trim_galore自动去除接头序列'
                }
            
            elif 'overrepresented sequences' in module_lower:
                issues['overrepresented'] = {
                    'severity': status,
                    'description': '过表达序列',
                    'solution': '检查是否存在污染或技术重复'
                }
            
            elif 'per tile sequence quality' in module_lower:
                issues['tile_quality'] = {
                    'severity': status,
                    'description': '测序芯片质量不均',
                    'solution': '可能是测序仪问题，通常不影响下游分析'
                }
            
            elif 'sequence length distribution' in module_lower:
                issues['length_distribution'] = {
                    'severity': status,
                    'description': '序列长度分布异常',
                    'solution': '检查文库构建和测序参数'
                }
    
    return issues

def generate_optimization_plan(issues):
    """生成质控优化方案"""
    print(f"\n=== 质控优化方案 ===")
    
    if not issues:
        print("✓ 数据质量良好，无需特殊优化")
        return
    
    # 按严重程度分类
    fail_issues = {k: v for k, v in issues.items() if v['severity'] == 'FAIL'}
    warn_issues = {k: v for k, v in issues.items() if v['severity'] == 'WARN'}
    
    print(f"\n严重问题 (FAIL): {len(fail_issues)} 个")
    for issue_name, issue_info in fail_issues.items():
        print(f"  ❌ {issue_info['description']}")
        print(f"     解决方案: {issue_info['solution']}")
    
    print(f"\n警告问题 (WARN): {len(warn_issues)} 个")
    for issue_name, issue_info in warn_issues.items():
        print(f"  ⚠️  {issue_info['description']}")
        print(f"     解决方案: {issue_info['solution']}")
    
    # 生成优化建议
    optimization_steps = []
    
    if 'adapter' in issues:
        optimization_steps.append("1. 使用trim_galore进行接头去除和质量修剪")
    
    if 'base_quality' in issues:
        optimization_steps.append("2. 设置质量阈值过滤低质量reads")
    
    if 'overrepresented' in issues:
        optimization_steps.append("3. 检查过表达序列是否为污染")
    
    if optimization_steps:
        print(f"\n推荐优化步骤:")
        for step in optimization_steps:
            print(f"  {step}")

def main():
    """主分析流程"""
    qc_dir = Path("analysis_results_1114/quality_control")
    
    if not qc_dir.exists():
        print(f"错误: 质控目录不存在: {qc_dir}")
        return
    
    html_files = list(qc_dir.glob("*.html"))
    
    if not html_files:
        print("未找到FastQC HTML报告")
        return
    
    all_issues = {}
    
    for html_file in html_files:
        modules = parse_fastqc_detailed(html_file)
        issues = analyze_qc_issues(modules)
        all_issues[html_file.name] = issues
    
    # 汇总所有样本的问题
    combined_issues = {}
    for sample_issues in all_issues.values():
        for issue_name, issue_info in sample_issues.items():
            if issue_name not in combined_issues:
                combined_issues[issue_name] = issue_info
            else:
                # 取更严重的状态
                if issue_info['severity'] == 'FAIL':
                    combined_issues[issue_name] = issue_info
    
    generate_optimization_plan(combined_issues)
    
    print(f"\n=== 分析完成 ===")
    print("\n下一步行动建议:")
    print("1. 根据上述优化方案改进质控流程")
    print("2. 重新运行质控分析验证改进效果")
    print("3. 确认质控通过后再进行下游分析")

if __name__ == "__main__":
    main()