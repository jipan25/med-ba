# 确保在正确的环境中
conda activate sc_analysis

# 获取当前 Conda 环境的根目录路径
CONDA_BASE=$(conda info --json | grep -o '"active_prefix":\s*"[^"]*' | cut -d'"' -f4)
MAKECONF="$CONDA_BASE/lib/R/etc/Makeconf"

# 备份原始配置文件 (以防万一)
cp "$MAKECONF" "$MAKECONF.bak.$(date +%s)"

# ----------------------------------------------------------------
# 关键步骤：重写 Makeconf
# 强制指定 C++ 标准，并强制包含 cairo/freetype 的头文件路径
# ----------------------------------------------------------------
cat > "$MAKECONF" <<EOF
# 使用 Conda 的编译器
CC = $CONDA_BASE/bin/x86_64-conda-linux-gnu-cc
CXX = $CONDA_BASE/bin/x86_64-conda-linux-gnu-c++
CXX11 = $CONDA_BASE/bin/x86_64-conda-linux-gnu-c++
CXX14 = $CONDA_BASE/bin/x86_64-conda-linux-gnu-c++
CXX17 = $CONDA_BASE/bin/x86_64-conda-linux-gnu-c++

# 强制所有 C++ 标准都使用 C++14 和 fPIC
CXXFLAGS = -std=gnu++14 -fPIC -g -O2
CXX11FLAGS = -std=gnu++14 -fPIC -g -O2
CXX14FLAGS = -std=gnu++14 -fPIC -g -O2
CXX17FLAGS = -std=gnu++14 -fPIC -g -O2

# 关键：将 Conda 的 include 目录加入预处理路径，解决 cairo-ft.h 找不到的问题
CPPFLAGS = -I$CONDA_BASE/include -I$CONDA_BASE/include/cairo -I$CONDA_BASE/include/freetype2 -DNDEBUG

# 链接器路径
LDFLAGS = -L$CONDA_BASE/lib -Wl,-rpath-link,$CONDA_BASE/lib
EOF

echo "✅ Makeconf 已强制覆写。编译路径已指向 Conda 环境。"
