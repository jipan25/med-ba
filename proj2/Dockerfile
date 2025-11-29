# *** 最终解决方案：使用 rocker/r-ver 基础镜像，然后手动安装 FULL Bioc 和拓展包 ***
# 使用 R 4.4.0 (对应 Bioc 3.22)，这个基础镜像体积小，拉取稳定。
FROM bioconductor/bioconductor_docker:latest

# 设置维护者和工作目录
LABEL maintainer="jipan25@gmail.com"
WORKDIR /home/rstudio

# 仅安装额外的系统依赖：主要是用于 R Markdown 报告生成的 LaTeX 工具
# Bioconductor 的其他系统依赖（如 libxml2-dev, libhdf5-dev）在这个基础镜像中已包含。
RUN apt-get update && apt-get install -y \
    # LaTeX 依赖：安装全套 TeXLive 包，以配合 tinytex::install_tinytex()
    texlive-full \
    # 常用工具（如果基础镜像没有，但通常是包含的）
    git \
    wget \
    # 清理以减小镜像
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# 复制 R 安装脚本到容器的临时目录
COPY install_packages.R /tmp/

# 执行安装脚本：安装 FULL Bioconductor + CRAN 拓展包，并安装 tinytex
# 警告：这一步耗时较长！
RUN Rscript /tmp/install_packages.R

# 清理：移除安装脚本
RUN rm /tmp/install_packages.R

# 设置容器启动时执行的默认命令 (例如启动 R 会话)
CMD ["R"]