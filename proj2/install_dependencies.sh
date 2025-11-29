#!/bin/bash

#定义虚拟环境名称

VENV_NAME=".venv"

echo "--- 正在创建和激活 Python 虚拟环境 ($VENV_NAME) ---"

# 1. 检查并创建虚拟环境

if [ ! -d "$VENV_NAME" ]; then
# 使用 python3 -m venv 创建虚拟环境
python3 -m venv $VENV_NAME
if [ $? -ne 0 ]; then
echo "错误：无法创建虚拟环境。请确保您安装了 python3-venv 包。"
exit 1
fi
echo "虚拟环境创建成功。"
else
echo "虚拟环境已存在。"
fi

 
 
PIP_PATH="./$VENV_NAME/bin/pip"

if [ ! -f "$PIP_PATH" ]; then
echo "错误：未找到虚拟环境中的 pip 可执行文件。请检查创建过程。"
exit 1
fi

echo "--- 正在虚拟环境中安装 Python 依赖：Scanpy 及其核心依赖 ---"
echo "--- 正在使用清华大学镜像源 (Tuna) 加速下载 ---"

# 使用虚拟环境中的 pip 安装所需的包，指定清华大学镜像源

$PIP_PATH install scanpy numpy pandas matplotlib \
-i https://pypi.tuna.tsinghua.edu.cn/simple

if [ $? -eq 0 ]; then
echo "--- 所有依赖安装成功！---"
echo "--- 接下来，如果您要在终端中运行 Python 脚本，请先手动激活环境："
echo "    source .venv/bin/activate"
echo "--- 当您完成工作后，输入 'deactivate' 即可退出虚拟环境。---"
else
echo "--- 依赖安装失败。请检查网络连接。---"
fi