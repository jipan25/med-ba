#!/bin/bash

# 1. -it: 交互式地运行容器，并分配一个伪终端
# 2. -v /home/ma/data:/data: 将您的本地数据目录挂载到容器内的 /data 目录。
# 3. -v /home/ma/ba:/ba: 将您的本地程序（脚本/代码）目录挂载到容器内的 /ba 目录。
# 4. bioconductor/bioconductor_docker:latest: 指定要使用的镜像名称。
# 5. /bin/bash: 在容器启动后执行 Bash shell，让您可以进入容器进行操作。

docker run -it \
    -v /home/ma/data:/data \
    -v /home/ma/ba:/ba \
    --cpus="12" \
    --name my_bio_container \
    -e R_BIOCPARALLEL_MAX_CORES=12 \
    -e R_AVAILABLE_CORES=12 \
    my-bio-project \
    /bin/bash

# 启动后，您将在容器内的命令行环境中。
# 此时：
# - 您的原始数据文件位于容器的 /data 目录下。
# - 您创建的 R 脚本文件位于容器的 /ba 目录下。