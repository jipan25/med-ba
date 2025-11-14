#!/bin/bash

# 服务器字体安装脚本
# 解决matplotlib中文字体问题

echo "=== 安装服务器字体支持 ==="

# 检查操作系统
if [[ "$OSTYPE" == "linux-gnu"* ]]; then
    echo "检测到Linux系统"
    
    # 检查是否支持中文
    if locale -a | grep -i zh > /dev/null; then
        echo "✅ 系统支持中文语言环境"
    else
        echo "⚠️  系统可能不支持中文，将安装字体包"
    fi
    
    # 安装字体包（Ubuntu/Debian）
    if command -v apt-get &> /dev/null; then
        echo "安装字体包..."
        sudo apt-get update
        sudo apt-get install -y fonts-dejavu fonts-liberation fontconfig
        
        # 安装中文字体包（可选）
        sudo apt-get install -y fonts-wqy-microhei fonts-wqy-zenhei
        
    # CentOS/RHEL
    elif command -v yum &> /dev/null; then
        echo "安装字体包..."
        sudo yum install -y dejavu-sans-fonts liberation-fonts
        sudo yum install -y wqy-microhei-fonts wqy-zenhei-fonts
    fi
    
    # 更新字体缓存
    echo "更新字体缓存..."
    sudo fc-cache -f -v
    
elif [[ "$OSTYPE" == "darwin"* ]]; then
    echo "检测到macOS系统，字体通常已安装"
else
    echo "未知系统类型: $OSTYPE"
fi

# 创建matplotlib字体配置文件
echo "配置matplotlib字体设置..."

# 创建matplotlib配置目录
mkdir -p ~/.config/matplotlib

# 创建字体配置文件
cat > ~/.config/matplotlib/matplotlibrc << 'EOF'
# 服务器字体配置
backend : Agg
font.family : sans-serif
font.sans-serif : DejaVu Sans, Liberation Sans, Arial, Helvetica, sans-serif
axes.unicode_minus : False

# 禁用交互式显示
interactive : False

# 图片设置
figure.dpi : 100
savefig.dpi : 300
savefig.bbox : tight
savefig.format : png
EOF

echo "✅ 字体配置完成"

# 测试字体配置
echo "测试字体配置..."
python3 -c "
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

# 创建测试图表
plt.figure(figsize=(8, 6))
x = np.linspace(0, 10, 100)
y = np.sin(x)
plt.plot(x, y)
plt.title('Font Test - English Title')
plt.xlabel('X Axis')
plt.ylabel('Y Axis')
plt.savefig('font_test.png')
print('✅ 字体测试图表已保存为 font_test.png')

# 检查可用字体
from matplotlib import font_manager
fonts = [f.name for f in font_manager.fontManager.ttflist if 'sans' in f.name.lower()]
print(f'可用无衬线字体: {set(fonts[:10])}')
"

echo ""
echo "=== 字体配置完成 ==="
echo "现在可以正常运行分析脚本:"
echo "python3 simplified_ba_analysis.py"
echo "python3 ba_analysis_pipeline.py"