#!/usr/bin/env python3
"""
胆道闭锁(BA)生物信息学分析流程 - Python版本
使用Python替代R进行差异表达分析
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from sklearn.decomposition import PCA
from statsmodels.stats.multitest import multipletests
import warnings
warnings.filterwarnings('ignore')

# 设置字体 - 服务器兼容版本
import matplotlib
matplotlib.use('Agg')  # 使用非交互式后端

# 尝试多种字体，优先使用服务器可用的字体
try:
    plt.rcParams['font.sans-serif'] = ['DejaVu Sans', 'Liberation Sans', 'Arial', 'Helvetica']
    plt.rcParams['axes.unicode_minus'] = False
    # 测试字体是否可用
    plt.figure()
    plt.title('Font Test')
    plt.close()
except:
    # 如果字体仍然有问题，使用默认设置
    plt.rcParams.update(plt.rcParamsDefault)
    print("警告: 使用默认字体设置")

def load_metadata(file_path):
    """加载样本元数据"""
    metadata = pd.read_excel(file_path)
    print("=== 样本元数据 ===")
    print(f"样本数量: {len(metadata)}")
    print(f"样本名称: {metadata['Sample name'].tolist()}")
    return metadata

def simulate_expression_data(metadata, n_genes=20000):
    """模拟基因表达数据"""
    np.random.seed(42)
    
    genes = [f'Gene_{i:05d}' for i in range(n_genes)]
    samples = metadata['Sample name'].tolist()
    
    # 创建表达矩阵
    expression_matrix = np.random.lognormal(mean=5, sigma=2, size=(n_genes, len(samples)))
    
    # 添加组间差异
    ba_indices = [i for i, sample in enumerate(samples) if 'RRV' in sample]
    nc_indices = [i for i, sample in enumerate(samples) if 'NC' in sample]
    
    # 在BA组中上调NETs相关基因
    nets_genes_indices = np.random.choice(n_genes, size=500, replace=False)
    expression_matrix[nets_genes_indices[:, None], ba_indices] *= 2.5
    
    # 在BA组中上调CD177相关基因
    cd177_genes_indices = np.random.choice(n_genes, size=300, replace=False)
    expression_matrix[cd177_genes_indices[:, None], ba_indices] *= 3.0
    
    expression_data = pd.DataFrame(expression_matrix, index=genes, columns=samples)
    
    print("=== 模拟表达数据 ===")
    print(f"基因数量: {n_genes}")
    print(f"数据形状: {expression_data.shape}")
    
    return expression_data

def differential_expression_analysis(expression_data):
    """差异表达分析 - Python实现"""
    print("\n=== 差异表达分析 ===")
    
    ba_samples = [col for col in expression_data.columns if 'RRV' in col]
    nc_samples = [col for col in expression_data.columns if 'NC' in col]
    
    results = []
    
    for gene in expression_data.index:
        ba_data = expression_data.loc[gene, ba_samples]
        nc_data = expression_data.loc[gene, nc_samples]
        
        # t检验
        t_stat, p_value = stats.ttest_ind(ba_data, nc_data)
        
        # 计算log2 fold change
        ba_mean = np.mean(ba_data)
        nc_mean = np.mean(nc_data)
        log2_fc = np.log2(ba_mean / nc_mean) if nc_mean > 0 else 0
        
        results.append({
            'gene': gene,
            'log2_fold_change': log2_fc,
            'p_value': p_value,
            'ba_mean': ba_mean,
            'nc_mean': nc_mean
        })
    
    de_results = pd.DataFrame(results)
    de_results['adj_p_value'] = multipletests(de_results['p_value'], method='fdr_bh')[1]
    
    # 筛选显著差异基因
    significant_genes = de_results[
        (de_results['adj_p_value'] < 0.05) & 
        (abs(de_results['log2_fold_change']) > 1)
    ]
    
    print(f"显著差异基因数量: {len(significant_genes)}")
    
    # 保存结果
    de_results.to_csv('differential_expression_results.csv', index=False)
    significant_genes.to_csv('significant_genes.csv', index=False)
    
    return de_results, significant_genes

def create_volcano_plot(de_results):
    """创建火山图"""
    plt.figure(figsize=(10, 8))
    
    # 设置颜色
    colors = []
    for _, row in de_results.iterrows():
        if row['adj_p_value'] < 0.05 and abs(row['log2_fold_change']) > 1:
            if row['log2_fold_change'] > 0:
                colors.append('red')  # 上调
            else:
                colors.append('blue')  # 下调
        else:
            colors.append('gray')  # 不显著
    
    plt.scatter(de_results['log2_fold_change'], 
                -np.log10(de_results['adj_p_value']), 
                c=colors, alpha=0.6, s=20)
    
    plt.axhline(y=-np.log10(0.05), color='red', linestyle='--', alpha=0.8)
    plt.axvline(x=1, color='red', linestyle='--', alpha=0.8)
    plt.axvline(x=-1, color='red', linestyle='--', alpha=0.8)
    
    plt.xlabel('log2 Fold Change')
    plt.ylabel('-log10(adj p-value)')
    plt.title('胆道闭锁 vs 正常对照 - 火山图')
    plt.grid(True, alpha=0.3)
    
    # 添加图例
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='red', label='上调基因'),
        Patch(facecolor='blue', label='下调基因'),
        Patch(facecolor='gray', label='不显著基因')
    ]
    plt.legend(handles=legend_elements)
    
    plt.tight_layout()
    plt.savefig('volcano_plot.png', dpi=300, bbox_inches='tight')
    plt.show()

def analyze_nets_genes(de_results):
    """分析NETs相关基因"""
    print("\n=== NETs相关基因分析 ===")
    
    # 模拟NETs相关基因（实际应根据文献定义）
    nets_genes = [f'Gene_{i:05d}' for i in np.random.choice(20000, 50, replace=False)]
    
    nets_results = de_results[de_results['gene'].isin(nets_genes)]
    
    if len(nets_results) > 0:
        print(f"找到 {len(nets_results)} 个NETs相关基因")
        
        # 统计
        upregulated = nets_results[nets_results['log2_fold_change'] > 1]
        downregulated = nets_results[nets_results['log2_fold_change'] < -1]
        
        print(f"上调基因: {len(upregulated)}")
        print(f"下调基因: {len(downregulated)}")
        print(f"平均log2FC: {nets_results['log2_fold_change'].mean():.2f}")
        
        # 保存结果
        nets_results.to_csv('nets_genes_analysis.csv', index=False)
    
    return nets_results

def main():
    """主函数"""
    print("开始胆道闭锁生物信息学分析...")
    
    # 加载元数据
    metadata = load_metadata('CRA007360.xlsx')
    
    # 模拟表达数据
    expression_data = simulate_expression_data(metadata)
    
    # 差异表达分析
    de_results, significant_genes = differential_expression_analysis(expression_data)
    
    # 创建火山图
    create_volcano_plot(de_results)
    
    # 分析NETs相关基因
    nets_results = analyze_nets_genes(de_results)
    
    print("\n=== 分析完成 ===")
    print("生成的文件:")
    print("- differential_expression_results.csv: 差异表达分析结果")
    print("- significant_genes.csv: 显著差异基因列表")
    print("- volcano_plot.png: 火山图")
    print("- nets_genes_analysis.csv: NETs相关基因分析")

if __name__ == "__main__":
    main()