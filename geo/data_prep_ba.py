import GEOparse
import pandas as pd
import numpy as np
from scipy import stats
from statsmodels.sandbox.stats.multicomp import multipletests
import os
import logging

# 配置日志
# 设置日志级别为 INFO，显示重要的进度信息
logging.basicConfig(level=logging.INFO, format='[%(levelname)s] %(message)s')

# --- 配置参数 ---
GEO_ID = 'GSE122340'
OUTPUT_DIR = 'geo_analysis_results'
GROUP_KEY = 'disease state' # 用于从样本元数据中识别分组信息的关键词
CONTROL_TERM = 'normal'
CASE_TERM = 'Biliary atresia'
LFC_THRESHOLD = 1.0  # log2FC 阈值
FDR_THRESHOLD = 0.05 # FDR 阈值

def find_expression_column(gsm_table):
    """
    尝试从GSM的表格数据中识别表达值列名。
    由于 GEOparse 默认解析可能失败，此函数通过尝试常见的列名来增强鲁棒性。
    """
    common_cols = ['Value', 'VALUE', 'Expression', 'Signal', 'V1', 'intensity', 'normalized_count']
    
    # 1. 检查常见名称
    for col in common_cols:
        if col in gsm_table.columns:
            logging.debug(f"找到常见的表达列名: {col}")
            return col

    # 2. 如果所有常见名称都不存在，尝试寻找第一个包含数字的列
    if not gsm_table.empty and len(gsm_table.columns) > 1:
        # 排除可能是ID的列（如'ID_REF'或第一个列）
        potential_data_cols = [c for c in gsm_table.columns if not c.upper().endswith('REF') and c != gsm_table.columns[0]]
        for col in potential_data_cols:
            try:
                # 检查前5行，看是否能转换为浮点数
                if gsm_table[col].head(5).astype(float).notna().all():
                    logging.warning(f"使用启发式方法找到的潜在表达列: {col}")
                    return col
            except (ValueError, TypeError):
                continue

    # 3. 如果所有方法都失败，则打印所有列名供调试
    logging.error(f"无法自动识别表达值列。GSM 表格的全部列名: {list(gsm_table.columns)}")
    raise ValueError("无法从 GEO 文件中提取有效的表达值列。请检查 DEBUG 日志中的列名。")

def load_and_process_data(geo_id):
    """加载GEO数据，处理表达矩阵并准备进行差异表达分析。"""
    print(f"\n--- 准备开始胆道闭锁生信复现之旅：GEO 数据的差异表达分析 ---")
    logging.info(f"--- 尝试加载 GEO 数据集: {geo_id} ---")

    # 确保输出目录存在，用于存放 GEOparse 下载的文件
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    try:
        # 1. 加载 GEO 数据
        # destdir 用于指定下载路径，方便管理
        gse = GEOparse.get_GEO(geo=geo_id, destdir=os.path.join(OUTPUT_DIR, geo_id), silent=True)
    except Exception as e:
        logging.error(f"加载 GEO 数据集失败。请检查网络连接或 GEO ID 是否正确。原始错误: {e}")
        # 修复 Bug: 确保返回的变量数量与接收变量数量一致 (data_matrix, group_df)
        return None, None

    logging.info(" -> 正在处理表达矩阵和样本分组...")
    
    data_matrix = None

    # 2. 识别表达值列并提取表达矩阵
    try:
        first_gsm_name = list(gse.gsms.keys())[0]
        first_gsm_table = gse.gsms[first_gsm_name].table
        
        # --- 鲁棒性检查: 如果 GSM 表格为空，则尝试直接读取 Series Matrix 文件 ---
        if first_gsm_table.empty or not list(first_gsm_table.columns):
            logging.warning("单个 GSM 表格内容为空 (列名为 [])。尝试从 Series Matrix 文件中直接加载数据。")
            
            matrix_path = os.path.join(OUTPUT_DIR, geo_id, f'{geo_id}_series_matrix.txt.gz')
            
            if not os.path.exists(matrix_path):
                 raise ValueError("GEOparse 无法解析 GSM 文件且 Series Matrix 文件未找到。")
                 
            # 找到数据起始行 (通常是包含 "ID_REF" 的行)
            data_start_line = 0
            with open(matrix_path, 'rt', encoding='utf-8') as f:
                for i, line in enumerate(f):
                    if line.startswith('"ID_REF"'):
                        data_start_line = i
                        break

            if data_start_line == 0:
                raise ValueError("无法在 Series Matrix 文件中找到数据起始行 (\"ID_REF\")。")

            logging.info(f" -> 发现数据起始于第 {data_start_line + 1} 行。")

            # 实际读取数据
            temp_matrix = pd.read_csv(
                matrix_path, 
                sep='\t', 
                skiprows=data_start_line, 
                index_col=0, 
                compression='gzip'
            )
            
            # 清理列名和索引名
            temp_matrix.columns = [col.strip('"') for col in temp_matrix.columns]
            temp_matrix.index = [idx.strip('"') for idx in temp_matrix.index]
            
            # 移除可能存在的尾部空列
            if temp_matrix.iloc[:, -1].isnull().all():
                temp_matrix = temp_matrix.iloc[:, :-1]
            
            # 确保数据是数字类型，并转置为 (Gene x Sample)
            data_matrix = temp_matrix.apply(pd.to_numeric, errors='coerce')
            data_matrix = data_matrix.T
            
            logging.info(f" -> 成功从 Series Matrix 文件中提取表达矩阵。形状: {data_matrix.shape}")
            
        else: # 如果 GSM 表格不为空，则使用原始的 GEOparse pivot 逻辑
            logging.debug(f"样本表格 ({first_gsm_name}) 中的所有列名:\n{list(first_gsm_table.columns)}")
            
            value_column = find_expression_column(first_gsm_table)
            logging.info(f" -> 识别到的表达值列名: {value_column}")

            # 3. 提取表达矩阵
            # pivot_samples 使用找到的列名来构建矩阵，转置后得到 (Gene x Sample)
            data_matrix = gse.pivot_samples(value_column).T
            logging.info(f" -> 成功提取表达矩阵。形状: {data_matrix.shape}")

    except ValueError as e:
        logging.error(f"--- 真实 GEO 数据加载失败 (ValueError: {e}) ---")
        print(f"--------------------------------------------------")
        logging.error(f"--- 程序因无法加载数据而中断 ---")
        return None, None
    except Exception as e:
        logging.error(f"提取表达矩阵时发生其他错误: {e}")
        return None, None

    # 4. 提取和清洗分组信息 (这部分逻辑保持不变，因为元数据通常是可靠的)
    groups = {}
    for gsm_name, gsm in gse.gsms.items():
        # GEOparse将characteristics解析为列表中的字符串
        characteristics = gsm.metadata.get('characteristics_ch1', [])
        
        # 查找包含 GROUP_KEY (disease state) 的特征
        group_info = next((char for char in characteristics if char.lower().startswith(GROUP_KEY)), None)
        
        if group_info:
            # 清理字符串，提取实际的值 (e.g., 'disease state: normal' -> 'normal')
            group_value = group_info.split(':')[-1].strip()
            groups[gsm_name] = group_value
            
    group_df = pd.Series(groups, name='Group')
    # 仅保留在表达矩阵中的样本
    group_df = group_df[group_df.index.isin(data_matrix.columns)] 
    
    # 检查分组是否有效
    if group_df.empty or len(group_df.unique()) < 2:
        logging.error("无法从元数据中提取有效的两组分组信息。")
        return None, None
    
    logging.info(f" -> 成功提取分组信息。共 {len(group_df)} 个样本，分组: {group_df.unique()}")

    # 确保表达矩阵的列和分组信息对齐
    data_matrix = data_matrix[group_df.index]
    
    return data_matrix, group_df

def perform_dea(data_matrix, group_df):
    """执行差异表达分析 (Differential Expression Analysis, DEA) 使用T-test。"""
    logging.info("--- 开始进行差异表达分析 (T-test) ---")

    # 分组样本
    control_samples = data_matrix.columns[group_df == CONTROL_TERM]
    case_samples = data_matrix.columns[group_df == CASE_TERM]

    logging.info(f" -> 控制组 ({CONTROL_TERM}) 样本数: {len(control_samples)}")
    logging.info(f" -> 疾病组 ({CASE_TERM}) 样本数: {len(case_samples)}")

    if len(control_samples) < 2 or len(case_samples) < 2:
        logging.error("分组样本数量不足 (每组至少需要2个样本) 无法进行T-test。")
        return None

    results = []

    # 遍历每个基因 (行)
    for gene in data_matrix.index:
        # 确保数据是浮点数类型
        case_data = data_matrix.loc[gene, case_samples].astype(float).dropna()
        control_data = data_matrix.loc[gene, control_samples].astype(float).dropna()
        
        # 检查数据是否足够
        if len(case_data) < 2 or len(control_data) < 2:
            t_stat, p_value, log2fc = np.nan, np.nan, np.nan
        else:
            # 独立样本 T-test (Welch's t-test，不等方差)
            t_stat, p_value = stats.ttest_ind(case_data, control_data, equal_var=False, nan_policy='omit')
            
            mean_case = case_data.mean()
            mean_control = control_data.mean()
            
            # 计算 Log2 Fold Change (LFC)
            # 为了避免 log(0) 错误，如果数据可能包含0，可以加上一个极小的数
            epsilon = 1e-5
            if mean_case >= 0 and mean_control >= 0:
                log2fc = np.log2((mean_case + epsilon) / (mean_control + epsilon))
            else:
                log2fc = np.nan

        results.append({
            'Gene': gene,
            'log2FC': log2fc,
            't_statistic': t_stat,
            'p_value': p_value
        })

    dea_results = pd.DataFrame(results).set_index('Gene')
    dea_results = dea_results.dropna(subset=['p_value', 'log2FC'])

    # 多重检验校正 (Benjamini-Hochberg False Discovery Rate)
    # 仅对非缺失的 p_value 进行校正
    p_values_to_correct = dea_results['p_value'].values
    reject, pvals_corrected, _, _ = multipletests(p_values_to_correct, method='fdr_bh')
    dea_results['FDR'] = pvals_corrected

    # 过滤和排序
    dea_results = dea_results.sort_values(by='FDR', ascending=True)

    logging.info(f"--- 差异表达分析完成。发现 {len(dea_results)} 个基因。 ---")
    
    return dea_results

def save_results(dea_results, geo_id):
    """保存差异表达基因结果。"""
    
    # 创建输出目录
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    # 筛选DEGs
    degs = dea_results[
        (dea_results['FDR'] < FDR_THRESHOLD) & 
        (np.abs(dea_results['log2FC']) >= LFC_THRESHOLD)
    ].copy() # 使用 .copy() 避免 SettingWithCopyWarning
    
    # 添加Up/Down列
    degs['Regulation'] = np.where(degs['log2FC'] > 0, 'Up', 'Down')
    
    output_path = os.path.join(OUTPUT_DIR, f'{geo_id}_DEA_full_results.xlsx')
    degs_path = os.path.join(OUTPUT_DIR, f'{geo_id}_DEGs_LFC{LFC_THRESHOLD}_FDR{FDR_THRESHOLD}.csv')

    # 保存完整结果 (Excel)
    dea_results.to_excel(output_path, index=True)
    logging.info(f" -> 完整DEA结果已保存到: {output_path}")

    # 保存筛选后的DEGs (CSV)
    degs[['log2FC', 'FDR', 'Regulation']].to_csv(degs_path, index=True)
    logging.info(f" -> {len(degs)} 个DEGs已保存到: {degs_path}")
    logging.info(f" -> 上调基因数: {len(degs[degs['Regulation'] == 'Up'])}")
    logging.info(f" -> 下调基因数: {len(degs[degs['Regulation'] == 'Down'])}")
    
    return degs

if __name__ == "__main__":
    # 运行主流程
    data_matrix, group_df = load_and_process_data(GEO_ID)

    if data_matrix is not None:
        dea_results = perform_dea(data_matrix, group_df)
        
        if dea_results is not None:
            degs = save_results(dea_results, GEO_ID)

            # 打印DEGs前几行作为摘要
            if not degs.empty:
                print("\n--- 差异表达基因 (Top 5) 摘要 ---")
                print(degs[['log2FC', 'FDR', 'Regulation']].head())
            else:
                print("\n--- 未发现满足阈值 (LFC>±1, FDR<0.05) 的差异表达基因 ---")
                
            print("\n--- 分析成功完成！请查看geo_analysis_results文件夹中的文件。 ---")