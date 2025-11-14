#!/usr/bin/env python3
"""
精确FastQC FAIL模块解析器
专门提取alt="[FAIL]"对应的具体模块名称
"""

import re
from pathlib import Path

def extract_fail_modules_precise(html_file):
    """精确提取FAIL模块名称"""
    print(f"\n=== 精确解析: {html_file.name} ===")
    
    with open(html_file, 'r', encoding='utf-8', errors='ignore') as f:
        content = f.read()
    
    # 查找FAIL标记及其上下文
    fail_patterns = [
        # 模式1: <td class="FAIL">模块名称</td>
        r'<td[^>]*class="FAIL"[^>]*>\s*([^<]+)\s*</td>',
        
        # 模式2: alt="[FAIL]"附近的模块名称
        r'alt="\[FAIL\]"[^>]*>\s*([^<]{10,100})\s*</td>',
        
        # 模式3: FAIL标记前后的模块名称
        r'([^>]{20,100})alt="\[FAIL\]"([^<]{20,100})',
        
        # 模式4: 表格行中的FAIL和模块名称
        r'<tr[^>]*>\s*<td[^>]*>\s*FAIL\s*</td>\s*<td[^>]*>([^<]+)</td>',
        
        # 模式5: 包含FAIL的完整表格行
        r'<tr[^>]*>\s*(?:<td[^>]*>[^<]*</td>\s*)*<td[^>]*>FAIL</td>\s*(?:<td[^>]*>[^<]*</td>\s*)*<td[^>]*>([^<]+)</td>',
    ]
    
    fail_modules = []
    
    for pattern in fail_patterns:
        matches = re.findall(pattern, content, re.IGNORECASE | re.DOTALL)
        for match in matches:
            if isinstance(match, tuple):
                # 处理多个捕获组
                for group in match:
                    if group and len(group.strip()) > 10:
                        module_name = group.strip()
                        # 清理HTML标签
                        module_name = re.sub(r'<[^>]+>', '', module_name)
                        if module_name and 'FAIL' not in module_name.upper():
                            fail_modules.append(module_name)
            else:
                # 单个捕获组
                if match and len(match.strip()) > 10:
                    module_name = match.strip()
                    module_name = re.sub(r'<[^>]+>', '', module_name)
                    if module_name and 'FAIL' not in module_name.upper():
                        fail_modules.append(module_name)
    
    # 去重和清理
    fail_modules = list(set([m.strip() for m in fail_modules if len(m.strip()) > 10]))
    
    # 如果还没有找到，尝试更宽松的模式
    if not fail_modules:
        # 查找包含"Per"或"Sequence"的文本块
        loose_pattern = r'([^>]{50,200}Per[^>]{50,200}alt="\[FAIL\]"[^<]{50,200})'
        loose_matches = re.findall(loose_pattern, content, re.IGNORECASE | re.DOTALL)
        
        for match in loose_matches:
            # 提取可能的模块名称
            module_candidates = re.findall(r'Per[^<>{}\[\]]{20,80}', match)
            for candidate in module_candidates:
                if len(candidate) > 15 and 'FAIL' not in candidate.upper():
                    fail_modules.append(candidate.strip())
    
    # 最终清理
    fail_modules = [m for m in fail_modules if len(m) > 10 and 'http' not in m and 'src=' not in m]
    
    fail_count = content.count('alt="[FAIL]"')
    print(f"发现FAIL标记数量: {fail_count}")
    print(f"识别FAIL模块: {len(fail_modules)} 个")
    
    if fail_modules:
        print("具体FAIL模块:")
        for i, module in enumerate(fail_modules, 1):
            print(f"  {i}. {module}")
    
    return fail_modules

def analyze_fail_impact(modules):
    """分析FAIL模块的影响"""
    print(f"\n=== FAIL影响分析 ===")
    
    critical_modules = [
        'Per base sequence quality',
        'Adapter Content',
        'Per sequence quality scores'
    ]
    
    warning_modules = [
        'Per base sequence content',
        'Per tile sequence quality',
        'Kmer Content'
    ]
    
    minor_modules = [
        'Sequence Length Distribution',
        'Overrepresented sequences',
        'Basic Statistics'
    ]
    
    critical_fails = []
    warning_fails = []
    minor_fails = []
    
    for module in modules:
        module_lower = module.lower()
        
        # 检查是否为关键模块
        is_critical = any(crit.lower() in module_lower for crit in critical_modules)
        is_warning = any(warn.lower() in module_lower for warn in warning_modules)
        is_minor = any(minor.lower() in module_lower for minor in minor_modules)
        
        if is_critical:
            critical_fails.append(module)
        elif is_warning:
            warning_fails.append(module)
        elif is_minor:
            minor_fails.append(module)
        else:
            # 未知模块，默认归为警告
            warning_fails.append(module)
    
    print(f"🔴 关键FAIL模块 ({len(critical_fails)} 个):")
    for module in critical_fails:
        print(f"   ❌ {module}")
    
    print(f"\n🟡 警告FAIL模块 ({len(warning_fails)} 个):")
    for module in warning_fails:
        print(f"   ⚠️  {module}")
    
    print(f"\n🟢 次要FAIL模块 ({len(minor_fails)} 个):")
    for module in minor_fails:
        print(f"   ℹ️  {module}")
    
    return {
        'critical': critical_fails,
        'warning': warning_fails,
        'minor': minor_fails
    }

def generate_optimization_recommendations(impact_analysis):
    """生成优化建议"""
    print(f"\n=== 质控优化建议 ===")
    
    recommendations = []
    
    # 关键FAIL的处理建议
    if impact_analysis['critical']:
        print("🔴 必须处理的严重问题:")
        
        for module in impact_analysis['critical']:
            if 'quality' in module.lower():
                recommendations.append({
                    'module': module,
                    'action': '质量修剪',
                    'tool': 'trim_galore',
                    'parameters': '--quality 20 --length 50',
                    'priority': '高'
                })
                print(f"   ❌ {module}: 使用trim_galore进行质量修剪")
            
            elif 'adapter' in module.lower():
                recommendations.append({
                    'module': module,
                    'action': '接头去除',
                    'tool': 'trim_galore',
                    'parameters': '--adapter AGATCGGAAGAGC',
                    'priority': '高'
                })
                print(f"   ❌ {module}: 必须去除接头序列")
    
    # 警告FAIL的处理建议
    if impact_analysis['warning']:
        print("\n🟡 建议处理的警告问题:")
        
        for module in impact_analysis['warning']:
            if 'content' in module.lower():
                recommendations.append({
                    'module': module,
                    'action': '碱基组成检查',
                    'tool': '人工检查',
                    'parameters': '通常可接受，但需注意',
                    'priority': '中'
                })
                print(f"   ⚠️  {module}: 检查碱基组成偏差")
    
    # 次要FAIL的处理建议
    if impact_analysis['minor']:
        print("\n🟢 可忽略的次要问题:")
        
        for module in impact_analysis['minor']:
            recommendations.append({
                'module': module,
                'action': '监控',
                'tool': '无需处理',
                'parameters': '通常不影响分析',
                'priority': '低'
            })
            print(f"   ℹ️  {module}: 通常可忽略")
    
    return recommendations

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
    
    print("=== FastQC FAIL模块精确分析 ===")
    
    all_recommendations = []
    
    for html_file in html_files:
        fail_modules = extract_fail_modules_precise(html_file)
        
        if fail_modules:
            impact = analyze_fail_impact(fail_modules)
            recommendations = generate_optimization_recommendations(impact)
            all_recommendations.extend(recommendations)
    
    # 生成总体优化方案
    print(f"\n=== 总体优化方案 ===")
    
    if not all_recommendations:
        print("✅ 未发现需要处理的FAIL模块")
        print("数据质量良好，可直接进行下游分析")
    else:
        # 按优先级排序
        high_priority = [r for r in all_recommendations if r['priority'] == '高']
        medium_priority = [r for r in all_recommendations if r['priority'] == '中']
        
        if high_priority:
            print("\n🔴 高优先级优化步骤:")
            for rec in high_priority:
                print(f"   1. {rec['action']}: {rec['tool']} {rec['parameters']}")
        
        if medium_priority:
            print("\n🟡 中优先级优化步骤:")
            for rec in medium_priority:
                print(f"   2. {rec['action']}: {rec['tool']}")
    
    print(f"\n=== 分析完成 ===")

if __name__ == "__main__":
    main()