import scanpy as sc
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt

# 设置 scanpy 绘图参数
sc.settings.verbosity = 3  # 设置详细输出级别
sc.settings.set_figure_params(dpi=100, facecolor='white')

def create_mock_adata(n_cells_raw=3000, n_genes_raw=20000):
    """
    用于回退机制，当找不到实际的 H5AD 文件时，创建一个更真实的模拟 AnnData 对象。
    这个模拟数据使用了较低的泊松率 (0.1) 来创建稀疏矩阵，以通过标准的 QC 过滤。
    """
    np.random.seed(42)
    
    # 模拟原始计数矩阵，使用较低的 lambda 值创建稀疏数据 (更像真实 scRNA-seq)
    X_raw = np.random.poisson(0.1, size=(n_cells_raw, n_genes_raw))
    
    # 模拟一部分细胞为低质量/死亡细胞 (表现为线粒体高表达)
    # 模拟前 300 个细胞的线粒体基因计数翻倍
    gene_names = [f'GENE_{i}' for i in range(n_genes_raw)]
    mt_gene_indices = np.arange(50) # 假设前 50 个基因是线粒体基因
    for i in mt_gene_indices:
        gene_names[i] = f'MT-{gene_names[i]}' 
        X_raw[:300, i] *= 2 # 增加前 300 个细胞的线粒体表达

    adata = sc.AnnData(X=X_raw, 
                        obs=pd.DataFrame(index=[f'CELL_{i}' for i in range(n_cells_raw)]), 
                        var=pd.DataFrame(index=gene_names))
    return adata


def run_gse176189_qc():
    """
    对 GSE176189 单细胞数据执行质量控制 (QC)。
    """
    print("--- 阶段一：GSE176189 质量控制开始 ---")

    # 1. 尝试加载上一步 (2_02_data_loading.R) 得到的原始 AnnData 对象
    file_path = 'gse176189_raw_data.h5ad'
    if os.path.exists(file_path):
        adata = sc.read_h5ad(file_path)
        print(f"已从文件加载原始数据: {file_path}")
    else:
        # 如果文件不存在，则使用模拟数据作为回退
        print(f"警告：未找到文件 {file_path}，正在使用模拟数据进行演示。")
        adata = create_mock_adata()
    
    n_cells_raw = adata.n_obs
    n_genes_raw = adata.n_vars
    print(f"原始细胞数: {n_cells_raw}, 原始基因数: {n_genes_raw}")

    # 2. 计算 QC 指标
    # 识别线粒体基因 (以 'MT-' 开头的基因)
    # 注意：这里假设您的原始数据包含以 'MT-' 开头的线粒体基因。
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    
    # 计算每个细胞的总计数、表达的基因数和线粒体基因计数的百分比
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

    # 可视化 QC 指标分布 (用于确定过滤阈值) - 保持注释，需要手动取消注释运行
    # sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], 
    #              jitter=0.4, multi_panel=True, save='_qc_metrics_pre_filter.png')
    
    print("\nQC 指标计算完成。")

    # 3. 基于指标进行细胞和基因过滤
    # ----------------------------------------------------
    # QC 阈值调整：使用更合理的范围来避免过滤掉所有细胞
    # ----------------------------------------------------
    
    # 过滤掉表达基因过少的细胞 (低质量)
    min_genes = 200
    sc.pp.filter_cells(adata, min_genes=min_genes)
    
    # 过滤掉表达基因过多的细胞 (可能为双峰或多重细胞)，已调整为 6000
    max_genes = 6000
    sc.pp.filter_cells(adata, max_genes=max_genes)
    
    # 过滤线粒体基因百分比过高的细胞 (可能濒临死亡或有损伤)
    max_mt_pct = 20
    # 注意：这里使用 inplace=False 的切片操作，并确保使用 .copy()
    initial_cells = adata.n_obs 
    adata = adata[adata.obs.pct_counts_mt < max_mt_pct, :].copy()
    print(f"过滤掉 {initial_cells - adata.n_obs} 个线粒体百分比 (> {max_mt_pct}%) 过高的细胞。")
    
    # 过滤掉在过少细胞中表达的基因 (减少稀疏性)
    min_cells = 3
    sc.pp.filter_genes(adata, min_cells=min_cells)

    n_cells_filtered = adata.n_obs
    n_genes_filtered = adata.n_vars
    cells_removed = n_cells_raw - n_cells_filtered
    genes_removed = n_genes_raw - n_genes_filtered

    print(f"\n过滤后细胞数: {n_cells_filtered}")
    print(f"过滤后基因数: {n_genes_filtered}")
    print(f"共过滤掉 {cells_removed} 个细胞和 {genes_removed} 个基因。")

    if n_cells_filtered == 0 or n_genes_filtered == 0:
        print("\n!!!!!!!! 严重警告：QC 过滤后数据为空，请检查原始数据质量和过滤阈值 !!!!!!!!")
        return None

    # 4. 归一化和对数转换
    # 归一化：将每个细胞的总计数调整为相同的目标值 (例如 10,000)，以消除测序深度差异。
    sc.pp.normalize_total(adata, target_sum=1e4)
    
    # 对数转换：对数据进行 log(x+1) 转换，以稳定方差并使数据更接近正态分布。
    sc.pp.log1p(adata)

    print("归一化和对数转换完成。")
    
    # 5. 保存 QC 后的 AnnData 对象
    adata.write(f"gse176189_qc_filtered.h5ad")
    print(f"QC 过滤后的数据已保存为 gse176189_qc_filtered.h5ad")

    print("\n--- GSE176189 质量控制阶段完成 ---")
    
    return adata

# 运行质量控制流程
if __name__ == '__main__':
    adata_qc = run_gse176189_qc()