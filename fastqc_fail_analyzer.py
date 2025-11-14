#!/usr/bin/env python3
"""
FastQC FAIL状态详细分析器
专门解析alt="[FAIL]"的质控问题
"""

import re
from pathlib import Path

def parse_fastqc_fail_status(html_file):
    """解析FastQC报告中的FAIL状态"""
    print(f"\n=== 详细解析FAIL状态: {html_file.name} ===")
    
    with open(html_file, 'r', encoding='utf-8', errors='ignore') as f:
        content = f.read()
    
    # 查找所有FAIL状态的图片标记
    fail_patterns = [
        r'alt="\[FAIL\]"',  # 标准FAIL标记
        r'alt="FAIL"',      # 可能的变体
        r'class="FAIL"',    # FAIL类
        r'>FAIL</',         # FAIL文本
    ]
    
    fail_matches = []
    for pattern in fail_patterns:
        matches = re.findall(pattern, content, re.IGNORECASE)
        if matches:
            fail_matches.extend(matches)
    
    # 查找FAIL对应的模块名称
    fail_modules = []
    
    # 查找模块名称的模式
    module_patterns = [
        r'alt="\[FAIL\]"[^>]*>\s*([^<]+)</td>',
        r'class="FAIL">([^<]+)</td>',
        r'>FAIL</td>\s*<td[^>]*>([^<]+)</td>',
    ]
    
    for pattern in module_patterns:
        matches = re.findall(pattern, content, re.IGNORECASE)
        if matches:
            fail_modules.extend([m.strip() for m in matches if len(m.strip()) > 5])
    
    # 如果直接查找失败，尝试查找FAIL附近的模块名称
    if not fail_modules and fail_matches:
        # 查找FAIL标记前后的文本内容
        context_pattern = r'([^>]{0,100})alt="\[FAIL\]"([^<]{0,100})'
        context_matches = re.findall(context_pattern, content, re.IGNORECASE)
        
        for before, after in context_matches:
            # 从前后文本中提取可能的模块名称
            module_candidates = re.findall(r'>([^<>]{10,50})</', before + after)
            for candidate in module_candidates:
                if len(candidate.strip()) > 5 and 'FAIL' not in candidate.upper():
                    fail_modules.append(candidate.strip())
    
    # 去重
    fail_modules = list(set(fail_modules))
    
    print(f"发现FAIL标记: {len(fail_matches)} 个")
    print(f"识别FAIL模块: {len(fail_modules)} 个")
    
    if fail_modules:
        print("具体FAIL模块:")
        for module in fail_modules:
            print(f"  ❌ {module}")
    
    return fail_modules

def analyze_fail_reasons(fail_modules, html_file):
    """分析FAIL原因并提供解决方案"""
    print(f"\n=== FAIL原因分析 ===")
    
    # 读取文件内容进行详细分析
    with open(html_file, 'r', encoding='utf-8', errors='ignore') as f:
        content = f.read()
    
    fail_analysis = {}
    
    for module in fail_modules:
        module_lower = module.lower()
        
        # 根据模块名称分析可能的原因
        if 'per base sequence quality' in module_lower:
            fail_analysis[module] = {
                'reason': '碱基质量分数过低或质量下降',
                'severity': '高',
                'impact': '影响比对准确性和基因定量',
                'solution': '使用trim_galore进行质量修剪，设置--quality 20',
                'parameters': '--quality 20 --length 50'
            }
        
        elif 'per sequence quality scores' in module_lower:
            fail_analysis[module] = {
                'reason': '序列质量分布异常',
                'severity': '中',
                'impact': '可能影响某些分析工具',
                'solution': '检查测序平台和文库构建质量',
                'parameters': '建议重新测序或使用质量控制工具'
            }
        
        elif 'per base sequence content' in module_lower:
            fail_analysis[module] = {
                'reason': '碱基组成存在系统性偏差',
                'severity': '中',
                'impact': '可能影响某些分析但通常可接受',
                'solution': '可能是测序起始偏差，通常可忽略',
                'parameters': '--clip_R1 10 --clip_R2 10 (去除前10bp)'
            }
        
        elif 'adapter content' in module_lower:
            fail_analysis[module] = {
                'reason': '检测到接头序列污染',
                'severity': '高',
                'impact': '严重影响比对和定量准确性',
                'solution': '必须去除接头序列',
                'parameters': '--adapter AGATCGGAAGAGC (Illumina通用接头)'
            }
        
        elif 'overrepresented sequences' in module_lower:
            fail_analysis[module] = {
                'reason': '存在过表达序列',
                'severity': '中',
                'impact': '可能表示污染或技术问题',
                'solution': '检查是否为rRNA污染或技术重复',
                'parameters': '使用kraken2进行物种鉴定'
            }
        
        elif 'per tile sequence quality' in module_lower:
            fail_analysis[module] = {
                'reason': '测序芯片质量不均',
                'severity': '低',
                'impact': '通常不影响下游分析',
                'solution': '测序仪问题，通常可接受',
                'parameters': '无需特殊处理'
            }
        
        elif 'sequence length distribution' in module_lower:
            fail_analysis[module] = {
                'reason': '序列长度分布异常',
                'severity': '中',
                'impact': '可能影响某些分析工具',
                'solution': '检查文库构建和测序参数',
                'parameters': '设置合适的长度过滤阈值'
            }
        
        else:
            # 未知模块的通用分析
            fail_analysis[module] = {
                'reason': '未知质控问题',
                'severity': '待定',
                'impact': '需要进一步分析',
                'solution': '查看详细FastQC报告',
                'parameters': '需要人工检查'
            }
    
    return fail_analysis

def generate_fail_optimization_plan(all_fail_analysis):
    """生成FAIL优化方案"""
    print(f"\n=== 质控FAIL优化方案 ===")
    
    if not all_fail_analysis:
        print("✓ 未发现FAIL模块，数据质量良好")
        return
    
    # 按严重程度分类
    high_severity = {}
    medium_severity = {}
    low_severity = {}
    
    for sample, analysis in all_fail_analysis.items():
        for module, info in analysis.items():
            if info['severity'] == '高':
                high_severity[f"{sample}: {module}"] = info
            elif info['severity'] == '中':
                medium_severity[f"{sample}: {module}"] = info
            else:
                low_severity[f"{sample}: {module}"] = info
    
    # 输出优化建议
    if high_severity:
        print(f"\n🔴 高严重性问题 ({len(high_severity)} 个):")
        for module_key, info in high_severity.items():
            print(f"\n❌ {module_key}")
            print(f"   原因: {info['reason']}")
            print(f"   影响: {info['impact']}")
            print(f"   解决方案: {info['solution']}")
            print(f"   参数建议: {info['parameters']}")
    
    if medium_severity:
        print(f"\n🟡 中等严重性问题 ({len(medium_severity)} 个):")
        for module_key, info in medium_severity.items():
            print(f"\n⚠️  {module_key}")
            print(f"   原因: {info['reason']}")
            print(f"   解决方案: {info['solution']}")
    
    if low_severity:
        print(f"\n🟢 低严重性问题 ({len(low_severity)} 个):")
        for module_key, info in low_severity.items():
            print(f"\nℹ️  {module_key}")
            print(f"   原因: {info['reason']}")
    
    # 生成具体的优化步骤
    print(f"\n=== 推荐优化步骤 ===")
    
    optimization_steps = []
    
    if any('adapter content' in key.lower() for key in high_severity.keys()):
        optimization_steps.append("1. 必须进行接头去除: trim_galore --adapter AGATCGGAAGAGC")
    
    if any('per base sequence quality' in key.lower() for key in high_severity.keys()):
        optimization_steps.append("2. 质量修剪: trim_galore --quality 20 --length 50")
    
    if any('overrepresented sequences' in key.lower() for key in medium_severity.keys()):
        optimization_steps.append("3. 检查过表达序列: 使用kraken2进行污染检测")
    
    if optimization_steps:
        for step in optimization_steps:
            print(f"  {step}")
    else:
        print("  无需特殊优化步骤")

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
    
    all_fail_analysis = {}
    
    for html_file in html_files:
        fail_modules = parse_fastqc_fail_status(html_file)
        
        if fail_modules:
            fail_analysis = analyze_fail_reasons(fail_modules, html_file)
            all_fail_analysis[html_file.name] = fail_analysis
    
    generate_fail_optimization_plan(all_fail_analysis)
    
    print(f"\n=== 分析完成 ===")
    print("\n下一步行动:")
    print("1. 根据FAIL分析结果优化质控流程")
    print("2. 实施推荐的优化步骤")
    print("3. 重新运行质控验证改进效果")

if __name__ == "__main__":
    main()