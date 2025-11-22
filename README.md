# med-ba
 ## qc质控
 1. trim_galore --paired -j 8 --quality 25 --length 75 --clip_R1 10 --clip_R2 10 --three_prime_clip_R1 10 --three_prime_clip_R2 10 --stringency 3 --trim-n --illumina --output_dir ~/soft-dev/data/results/optimized_trimmed/ CRR524834_f1.fq.gz CRR524834_r2.fq.gz

 2. trim_galore --paired -j 8   --quality 25   --length 75   --clip_R1 10   --clip_R2 10   --three_prime_clip_R1 10   --three_prime_clip_R2 10   --stringency 3   --trim-n   --illumina   --output_dir ~/soft-dev/data/results/final_trimmed/   CRR524835_f1.fq.gz CRR524835_r2.fq.gz

 3. - fastqc -t 8 -o ~/soft-dev/data/results/fastqc/ ~/soft-dev/data/results/final_trimmed/*.fq.gz
  - fastqc -t 8 -o ~/soft-dev/data/results/fastqc/ ~/soft-dev/data/results/optimized_trimmed/*.fq.gz 
 4. multiqc -i ~/soft-dev/data/results/fastqc/ ~/soft-dev/data/results/fastqc/ -o ~/soft-dev/data/results/multiqc/
 5. 

 ## 下载小鼠转录组（mm10）和注释文件
##### 创建参考转录组目录
mkdir -p ~/soft-dev/data/reference/mouse/

##### 下载小鼠转录组（mm10）和注释文件
cd ~/soft-dev/data/reference/mouse/

#####  下载转录组序列（不是基因组！）
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.transcripts.fa.gz
gunzip gencode.vM25.transcripts.fa.gz

##### 下载GTF注释文件
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.annotation.gtf.gz
gunzip gencode.vM25.annotation.gtf.gz

2. salmon index -t gencode.vM25.transcripts.fa \
  -i mm10_salmon_index \
  --gencode -p 8

* grep -v "^#" gencode.vM25.annotation.gtf | \
  awk '$3 == "transcript" {split($10, gene, "\""); split($12, tx, "\""); print tx[2], gene[2]}' > \
  transcript_to_gene_map.txt 


## 细胞定量分析

#### 为每个样本运行alevin定量
mkdir -p ~/soft-dev/data/results/alevin/

#### 处理CRR524834 (NC_5d对照组)
salmon alevin -l ISR \
  -i ~/soft-dev/data/reference/mouse/mm10_salmon_index \
  -1 ~/soft-dev/data/results/final_trimmed/CRR524834_f1_val_1.fq.gz \
  -2 ~/soft-dev/data/results/final_trimmed/CRR524834_r2_val_2.fq.gz \
  -o ~/soft-dev/data/results/alevin/CRR524834 \
  -p 12 \
  --chromiumV3 \
  --tgMap ~/soft-dev/data/reference/mouse/transcript_to_gene_map.txt

#### 处理CRR524835 (RRV_5d实验组)  
salmon alevin -l ISR \
  -i ~/soft-dev/data/reference/mouse/mm10_salmon_index \
  -1 ~/soft-dev/data/results/final_trimmed/CRR524835_f1_val_1.fq.gz \
  -2 ~/soft-dev/data/results/final_trimmed/CRR524835_r2_val_2.fq.gz \
  -o ~/soft-dev/data/results/alevin/CRR524835 \
  -p 12 \
  --chromiumV3 \
  --tgMap ~/soft-dev/data/reference/mouse/transcript_to_gene_map.txt


## gemini 3 告知不需要trim_galore ，直接运行
### 1. 创建结果目录
mkdir -p ~/soft-dev/data/results/alevin_raw/CRR524834
(strict_qc_env) ma@ma-MS:~/soft-dev/data/test_chem/res_v3$ salmon alevin -l ISR \
  -i ~/soft-dev/data/reference/mouse/mm10_salmon_index \
  -1 ~/soft-dev/data/med-data/CRA007360/CRR524834/CRR524834_f1.fq.gz \
  -2 ~/soft-dev/data/med-data/CRA007360/CRR524834/CRR524834_r2.fq.gz \
  -o ~/soft-dev/data/results/alevin_raw/CRR524834 \
  -p 12 \
  --chromiumV3 \
  --tgMap ~/soft-dev/data/reference/mouse/transcript_to_gene_map.txt \
  --dumpFeatures \
  --expectCells 5000

salmon alevin -l ISR \
  -i ~/soft-dev/data/reference/mouse/mm10_salmon_index \
  -1 ~/soft-dev/data/med-data/CRA007360/CRR524835/CRR524835_f1.fq.gz \
  -2 ~/soft-dev/data/med-data/CRA007360/CRR524835/CRR524835_r2.fq.gz \
  -o ~/soft-dev/data/results/alevin_raw/CRR524835 \
  -p 12 \
  --chromiumV3 \
  --tgMap ~/soft-dev/data/reference/mouse/transcript_to_gene_map.txt \
  --dumpFeatures \
  --expectCells 5000

### 2. 结论：数据处理成功！
高比对率确认：两个样本的比对率都稳定在 67% - 68% 左右。这不仅远高于我们最初的 0.07%，也符合高质量单细胞 RNA-seq 数据的标准（通常 40%-75% 都是可接受的）。

方法确认：这彻底证明了：

trim_galore 是错误的。

直接使用原始数据 (Raw data) + 正确的 chemistry flag (V3，尽管我们推断是 V2) 是正确的策略。

### 3. # 文件名: 01_data_import_and_qc.R
Rscript 01_data_import_and_qc.R
 
 由于原来  conda 报错 重新初始化R 语言环境
 切换到 conda create -n sc_analysis r=4.3 -c conda-forge

conda activate sc_analysis