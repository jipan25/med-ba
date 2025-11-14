#!/bin/bash

# 胆道闭锁分析项目环境配置脚本
# 适用于服务器环境

echo "=== 胆道闭锁生物信息学分析环境配置 ==="

# 检查Python版本
python_version=$(python3 -c "import sys; print('.'.join(map(str, sys.version_info[:2])))" 2>/dev/null)
if [ $? -ne 0 ]; then
    echo "错误: 未找到Python3，请先安装Python 3.8+"
    exit 1
fi

echo "检测到Python版本: $python_version"

# 检查pip是否可用
if ! command -v pip3 &> /dev/null; then
    echo "错误: 未找到pip3，请先安装pip"
    exit 1
fi

echo "=== 安装Python依赖包 ==="

# 升级pip
pip3 install --upgrade pip

# 安装核心依赖
pip3 install pandas numpy scipy statsmodels scikit-learn matplotlib seaborn openpyxl xlrd

# 检查安装是否成功
if python3 -c "import pandas, numpy, scipy, matplotlib, seaborn; print('依赖包安装成功')" 2>/dev/null; then
    echo "✅ 所有依赖包安装成功"
else
    echo "❌ 部分依赖包安装失败，请检查错误信息"
    exit 1
fi

echo "=== 环境验证 ==="

# 验证关键功能
python3 -c "
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from sklearn.decomposition import PCA
from statsmodels.stats.multitest import multipletests

print('✅ 核心包导入成功')

# 测试基本功能
data = pd.DataFrame({'A': [1,2,3], 'B': [4,5,6]})
print('✅ pandas功能正常')

arr = np.array([1,2,3])
print('✅ numpy功能正常')

# 测试统计功能
test_result = stats.ttest_ind([1,2,3], [4,5,6])
print('✅ scipy统计功能正常')

print('🎉 环境配置完成！可以运行分析脚本了。')
"

echo ""
echo "=== 使用说明 ==="
echo "1. 运行完整分析: python3 ba_analysis_pipeline.py"
echo "2. 运行简化分析: python3 simplified_ba_analysis.py"
echo "3. 运行Python版本: python3 python_ba_analysis.py"
echo ""
echo "分析结果将保存在当前目录下"