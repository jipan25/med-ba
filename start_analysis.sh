#!/bin/bash

# 胆道闭锁分析启动脚本
# 使用现有的 ba_analysis 环境

echo "=== 启动胆道闭锁真实数据分析 ==="

# 检查环境是否激活
if [[ "$CONDA_DEFAULT_ENV" != "ba_analysis" ]]; then
    echo "激活 ba_analysis 环境..."
    conda activate ba_analysis
fi

# 检查Python版本
echo "Python版本: $(python --version)"
echo "当前环境: $CONDA_DEFAULT_ENV"

# 安装必要的包
echo "安装分析依赖包..."
pip install -r requirements_ba_analysis.txt

# 运行分析
echo "开始数据分析..."
python real_data_analysis_pipeline.py

echo "=== 分析完成 ==="
echo "结果保存在: analysis_results/ 目录"