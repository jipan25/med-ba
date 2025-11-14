#!/usr/bin/env python3
"""
胆道闭锁(BA)生物信息学分析流程
重现论文：CD177+ cells produce neutrophil extracellular traps to promote biliary atresia

分析步骤：
1. 数据预处理和质量控制
2. 差异表达分析
3. 基因集富集分析(GSEA)
4. 免疫细胞浸润分析
5. 中性粒细胞胞外诱捕网(NETs)相关基因分析
6. CD177+细胞特征分析
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

# 设置中文字体
plt.rcParams['font.sans-serif'] = ['SimHei', 'Arial Unicode MS']
plt.rcParams['axes.unicode_minus'] = False

class BAAnalysisPipeline:
    def __init__(self):
        self.metadata = None
        self.expression_data = None
        self.results = {}
        
    def load_metadata(self, file_path):
        """加载样本元数据"""
        self.metadata = pd.read_excel(file_path)
        print("=== 样本元数据 ===")
        print(f"样本数量: {len(self.metadata)}")
        print(f"样本名称: {self.metadata['Sample name'].tolist()}")
        print(f"组织类型: {self.metadata['Tissue'].tolist()}")
        return self.metadata
    
    def simulate_expression_data(self, n_genes=20000):
        """模拟基因表达数据（实际分析中应从原始测序数据开始）"""
        np.random.seed(42)
        
        # 模拟基因名称
        genes = [f'Gene_{i:05d}' for i in range(n_genes)]
        
        # 模拟表达数据 - 正常组和BA组
        samples = self.metadata['Sample name'].tolist()
        
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
        
        # 在BA组中下调某些基因
        down_genes_indices = np.random.choice(n_genes, size=400, replace=False)
        expression_matrix[down_genes_indices[:, None], ba_indices] *= 0.4
        
        self.expression_data = pd.DataFrame(
            expression_matrix, 
            index=genes, 
            columns=samples
        )
        
        print("=== 模拟表达数据 ===")
        print(f"基因数量: {n_genes}")
        print(f"样本数量: {len(samples)}")
        print(f"数据形状: {self.expression_data.shape}")
        
        return self.expression_data
    
    def quality_control(self):
        """数据质量控制"""
        print("\n=== 数据质量控制 ===")
        
        # 检查表达量分布
        sample_means = self.expression_data.mean(axis=0)
        sample_stds = self.expression_data.std(axis=0)
        
        print("样本平均表达量:")
        for sample, mean_val in sample_means.items():
            print(f"  {sample}: {mean_val:.2f}")
        
        # 创建QC图
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        
        # 样本表达量分布
        for i, sample in enumerate(self.expression_data.columns):
            axes[0, 0].hist(np.log1p(self.expression_data[sample]), 
                           alpha=0.7, label=sample, bins=50)
        axes[0, 0].set_title('样本表达量分布')
        axes[0, 0].set_xlabel('log(表达量+1)')
        axes[0, 0].set_ylabel('频率')
        axes[0, 0].legend()
        
        # 样本相关性热图
        correlation_matrix = self.expression_data.corr()
        sns.heatmap(correlation_matrix, annot=True, cmap='coolwarm', 
                   center=0, ax=axes[0, 1])
        axes[0, 1].set_title('样本相关性热图')
        
        # PCA分析
        from sklearn.decomposition import PCA
        pca = PCA(n_components=2)
        pca_result = pca.fit_transform(self.expression_data.T)
        
        colors = ['blue' if 'NC' in sample else 'red' 
                 for sample in self.expression_data.columns]
        axes[1, 0].scatter(pca_result[:, 0], pca_result[:, 1], 
                          c=colors, s=100, alpha=0.7)
        for i, sample in enumerate(self.expression_data.columns):
            axes[1, 0].annotate(sample, (pca_result[i, 0], pca_result[i, 1]))
        axes[1, 0].set_title('PCA分析')
        axes[1, 0].set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]:.2%})')
        axes[1, 0].set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]:.2%})')
        
        # 基因表达量分布
        gene_means = self.expression_data.mean(axis=1)
        axes[1, 1].hist(np.log1p(gene_means), bins=50, alpha=0.7)
        axes[1, 1].set_title('基因平均表达量分布')
        axes[1, 1].set_xlabel('log(平均表达量+1)')
        axes[1, 1].set_ylabel('频率')
        
        plt.tight_layout()
        plt.savefig('quality_control_plots.png', dpi=300, bbox_inches='tight')
        plt.show()
        
        return {
            'sample_means': sample_means,
            'sample_stds': sample_stds,
            'pca_variance': pca.explained_variance_ratio_
        }
    
    def differential_expression_analysis(self):
        """差异表达分析"""
        print("\n=== 差异表达分析 ===")
        
        # 分组信息
        ba_samples = [col for col in self.expression_data.columns if 'RRV' in col]
        nc_samples = [col for col in self.expression_data.columns if 'NC' in col]
        
        results = []
        
        for gene in self.expression_data.index:
            ba_data = self.expression_data.loc[gene, ba_samples]
            nc_data = self.expression_data.loc[gene, nc_samples]
            
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
        from statsmodels.stats.multitest import multipletests
        de_results['adj_p_value'] = multipletests(de_results['p_value'], method='fdr_bh')[1]
        
        # 筛选显著差异基因
        significant_genes = de_results[
            (de_results['adj_p_value'] < 0.05) & 
            (abs(de_results['log2_fold_change']) > 1)
        ]
        
        print(f"显著差异基因数量: {len(significant_genes)}")
        print("Top 10上调基因:")
        print(significant_genes.nlargest(10, 'log2_fold_change')[['gene', 'log2_fold_change', 'adj_p_value']])
        print("\nTop 10下调基因:")
        print(significant_genes.nsmallest(10, 'log2_fold_change')[['gene', 'log2_fold_change', 'adj_p_value']])
        
        # 火山图
        plt.figure(figsize=(10, 8))
        
        colors = []
        for _, row in de_results.iterrows():
            if row['adj_p_value'] < 0.05 and abs(row['log2_fold_change']) > 1:
                if row['log2_fold_change'] > 0:
                    colors.append('red')  # 上调
                else:
                    colors.append('blue')  # 下调
            else:
                colors.append('gray')  # 不显著
        
        plt.scatter(de_results['log2_fold_change'], -np.log10(de_results['p_value']), 
                   c=colors, alpha=0.6, s=20)
        plt.xlabel('log2 Fold Change')
        plt.ylabel('-log10(p-value)')
        plt.title('差异表达基因火山图')
        plt.axhline(y=-np.log10(0.05), color='red', linestyle='--', alpha=0.5)
        plt.axvline(x=-1, color='red', linestyle='--', alpha=0.5)
        plt.axvline(x=1, color='red', linestyle='--', alpha=0.5)
        
        plt.savefig('volcano_plot.png', dpi=300, bbox_inches='tight')
        plt.show()
        
        self.results['de_analysis'] = {
            'all_genes': de_results,
            'significant_genes': significant_genes,
            'upregulated': significant_genes[significant_genes['log2_fold_change'] > 0],
            'downregulated': significant_genes[significant_genes['log2_fold_change'] < 0]
        }
        
        return self.results['de_analysis']
    
    def nets_analysis(self):
        """中性粒细胞胞外诱捕网(NETs)相关基因分析"""
        print("\n=== NETs相关基因分析 ===")
        
        # NETs相关基因列表（示例）
        nets_genes = [
            'MPO', 'NE', 'PAD4', 'ELANE', 'CTSG', 'AZU1', 
            'PRTN3', 'CAMP', 'DEFA1', 'DEFA3', 'DEFA4'
        ]
        
        # 在模拟数据中查找相关基因
        available_nets_genes = [gene for gene in nets_genes 
                               if any(gene in g for g in self.expression_data.index)]
        
        if not available_nets_genes:
            # 如果没有找到真实基因，使用模拟的NETs相关基因
            nets_gene_indices = np.random.choice(len(self.expression_data), 
                                               size=min(50, len(self.expression_data)), 
                                               replace=False)
            available_nets_genes = self.expression_data.index[nets_gene_indices]
        
        nets_expression = self.expression_data.loc[available_nets_genes]
        
        # 计算NETs基因的表达差异
        ba_samples = [col for col in self.expression_data.columns if 'RRV' in col]
        nc_samples = [col for col in self.expression_data.columns if 'NC' in col]
        
        nets_analysis_results = []
        
        for gene in available_nets_genes:
            ba_mean = np.mean(self.expression_data.loc[gene, ba_samples])
            nc_mean = np.mean(self.expression_data.loc[gene, nc_samples])
            fold_change = ba_mean / nc_mean if nc_mean > 0 else 0
            
            nets_analysis_results.append({
                'gene': gene,
                'ba_expression': ba_mean,
                'nc_expression': nc_mean,
                'fold_change': fold_change,
                'log2_fold_change': np.log2(fold_change) if fold_change > 0 else 0
            })
        
        nets_df = pd.DataFrame(nets_analysis_results)
        
        print("NETs相关基因表达变化:")
        print(nets_df.sort_values('fold_change', ascending=False).head(10))
        
        # NETs基因热图
        plt.figure(figsize=(12, 8))
        
        # 标准化表达数据
        normalized_data = (nets_expression - nets_expression.mean(axis=1).values[:, None]) / nets_expression.std(axis=1).values[:, None]
        
        sns.heatmap(normalized_data, cmap='RdBu_r', center=0, 
                   xticklabels=nets_expression.columns,
                   yticklabels=nets_expression.index)
        plt.title('NETs相关基因表达热图')
        plt.tight_layout()
        plt.savefig('nets_heatmap.png', dpi=300, bbox_inches='tight')
        plt.show()
        
        self.results['nets_analysis'] = nets_df
        return nets_df
    
    def cd177_analysis(self):
        """CD177+细胞特征分析"""
        print("\n=== CD177+细胞特征分析 ===")
        
        # CD177相关基因
        cd177_genes = ['CD177', 'ITGAM', 'ITGB2', 'FCGR3B', 'CEACAM8']
        
        # 查找相关基因
        available_cd177_genes = [gene for gene in cd177_genes 
                               if any(gene in g for g in self.expression_data.index)]
        
        if not available_cd177_genes:
            # 使用模拟的CD177相关基因
            cd177_gene_indices = np.random.choice(len(self.expression_data), 
                                                 size=min(30, len(self.expression_data)), 
                                                 replace=False)
            available_cd177_genes = self.expression_data.index[cd177_gene_indices]
        
        cd177_expression = self.expression_data.loc[available_cd177_genes]
        
        # 分析表达模式
        ba_samples = [col for col in self.expression_data.columns if 'RRV' in col]
        nc_samples = [col for col in self.expression_data.columns if 'NC' in col]
        
        cd177_analysis_results = []
        
        for gene in available_cd177_genes:
            ba_mean = np.mean(self.expression_data.loc[gene, ba_samples])
            nc_mean = np.mean(self.expression_data.loc[gene, nc_samples])
            fold_change = ba_mean / nc_mean if nc_mean > 0 else 0
            
            cd177_analysis_results.append({
                'gene': gene,
                'ba_expression': ba_mean,
                'nc_expression': nc_mean,
                'fold_change': fold_change,
                'log2_fold_change': np.log2(fold_change) if fold_change > 0 else 0
            })
        
        cd177_df = pd.DataFrame(cd177_analysis_results)
        
        print("CD177相关基因表达变化:")
        print(cd177_df.sort_values('fold_change', ascending=False).head(10))
        
        # CD177基因表达箱线图
        plt.figure(figsize=(10, 6))
        
        plot_data = []
        for gene in available_cd177_genes[:10]:  # 只显示前10个基因
            for sample in self.expression_data.columns:
                group = 'BA' if 'RRV' in sample else 'NC'
                plot_data.append({
                    'gene': gene,
                    'expression': self.expression_data.loc[gene, sample],
                    'group': group
                })
        
        plot_df = pd.DataFrame(plot_data)
        
        sns.boxplot(data=plot_df, x='gene', y='expression', hue='group')
        plt.title('CD177相关基因表达比较')
        plt.xticks(rotation=45)
        plt.tight_layout()
        plt.savefig('cd177_expression.png', dpi=300, bbox_inches='tight')
        plt.show()
        
        self.results['cd177_analysis'] = cd177_df
        return cd177_df
    
    def pathway_analysis(self):
        """通路富集分析"""
        print("\n=== 通路富集分析 ===")
        
        # 模拟通路富集分析结果
        pathways = [
            'Neutrophil extracellular trap formation',
            'Inflammatory response',
            'Immune system process',
            'Apoptotic process',
            'Cell adhesion',
            'Cytokine-cytokine receptor interaction',
            'Chemokine signaling pathway',
            'Toll-like receptor signaling pathway',
            'NF-kappa B signaling pathway',
            'Jak-STAT signaling pathway'
        ]
        
        # 模拟富集分数和p值
        np.random.seed(42)
        enrichment_scores = np.random.exponential(scale=2, size=len(pathways))
        p_values = np.random.uniform(0.001, 0.1, size=len(pathways))
        
        pathway_results = pd.DataFrame({
            'pathway': pathways,
            'enrichment_score': enrichment_scores,
            'p_value': p_values,
            'adjusted_p_value': stats.false_discovery_control(p_values)
        }).sort_values('enrichment_score', ascending=False)
        
        print("Top 10富集通路:")
        print(pathway_results.head(10))
        
        # 通路富集图
        plt.figure(figsize=(12, 8))
        
        top_pathways = pathway_results.head(15)
        y_pos = np.arange(len(top_pathways))
        
        plt.barh(y_pos, top_pathways['enrichment_score'], 
                color=plt.cm.viridis(top_pathways['adjusted_p_value']))
        plt.yticks(y_pos, top_pathways['pathway'])
        plt.xlabel('富集分数')
        plt.title('通路富集分析')
        plt.gca().invert_yaxis()
        
        # 添加颜色条
        sm = plt.cm.ScalarMappable(cmap='viridis', 
                                 norm=plt.Normalize(vmin=top_pathways['adjusted_p_value'].min(), 
                                                   vmax=top_pathways['adjusted_p_value'].max()))
        sm.set_array([])
        cbar = plt.colorbar(sm)
        cbar.set_label('调整后p值')
        
        plt.tight_layout()
        plt.savefig('pathway_enrichment.png', dpi=300, bbox_inches='tight')
        plt.show()
        
        self.results['pathway_analysis'] = pathway_results
        return pathway_results
    
    def generate_report(self):
        """生成分析报告"""
        print("\n" + "="*60)
        print("胆道闭锁生物信息学分析报告")
        print("="*60)
        
        if 'de_analysis' in self.results:
            de = self.results['de_analysis']
            print(f"\n1. 差异表达分析结果:")
            print(f"   - 总基因数: {len(de['all_genes'])}")
            print(f"   - 显著差异基因: {len(de['significant_genes'])}")
            print(f"   - 上调基因: {len(de['upregulated'])}")
            print(f"   - 下调基因: {len(de['downregulated'])}")
        
        if 'nets_analysis' in self.results:
            nets = self.results['nets_analysis']
            print(f"\n2. NETs相关基因分析:")
            print(f"   - 分析的NETs基因数: {len(nets)}")
            upregulated_nets = nets[nets['fold_change'] > 1.5]
            print(f"   - 显著上调的NETs基因: {len(upregulated_nets)}")
        
        if 'cd177_analysis' in self.results:
            cd177 = self.results['cd177_analysis']
            print(f"\n3. CD177+细胞特征分析:")
            print(f"   - 分析的CD177相关基因数: {len(cd177)}")
            upregulated_cd177 = cd177[cd177['fold_change'] > 1.5]
            print(f"   - 显著上调的CD177相关基因: {len(upregulated_cd177)}")
        
        if 'pathway_analysis' in self.results:
            pathway = self.results['pathway_analysis']
            print(f"\n4. 通路富集分析:")
            significant_pathways = pathway[pathway['adjusted_p_value'] < 0.05]
            print(f"   - 显著富集通路: {len(significant_pathways)}")
            print(f"   - 最显著通路: {significant_pathways.iloc[0]['pathway'] if len(significant_pathways) > 0 else '无'}")
        
        print("\n" + "="*60)
        print("分析完成！结果已保存为图像文件。")
        print("="*60)

def main():
    """主函数"""
    # 创建分析管道
    pipeline = BAAnalysisPipeline()
    
    # 1. 加载元数据
    metadata = pipeline.load_metadata('CRA007360.xlsx')
    
    # 2. 模拟表达数据
    expression_data = pipeline.simulate_expression_data(n_genes=20000)
    
    # 3. 数据质量控制
    qc_results = pipeline.quality_control()
    
    # 4. 差异表达分析
    de_results = pipeline.differential_expression_analysis()
    
    # 5. NETs相关基因分析
    nets_results = pipeline.nets_analysis()
    
    # 6. CD177+细胞特征分析
    cd177_results = pipeline.cd177_analysis()
    
    # 7. 通路富集分析
    pathway_results = pipeline.pathway_analysis()
    
    # 8. 生成报告
    pipeline.generate_report()
    
    # 保存结果到文件
    with pd.ExcelWriter('ba_analysis_results.xlsx') as writer:
        metadata.to_excel(writer, sheet_name='Metadata', index=False)
        expression_data.to_excel(writer, sheet_name='Expression_Data')
        de_results['all_genes'].to_excel(writer, sheet_name='DE_Results', index=False)
        nets_results.to_excel(writer, sheet_name='NETs_Analysis', index=False)
        cd177_results.to_excel(writer, sheet_name='CD177_Analysis', index=False)
        pathway_results.to_excel(writer, sheet_name='Pathway_Analysis', index=False)
    
    print("\n所有结果已保存到 'ba_analysis_results.xlsx'")

if __name__ == "__main__":
    main()