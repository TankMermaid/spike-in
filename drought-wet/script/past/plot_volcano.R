#!/usr/bin/env Rscript
# 
# Copyright 2016-2018 Yong-Xin Liu <metagenome@126.com>

# If used this script, please cited:
#   Zhang, J., Zhang, N., Liu, Y.X., Zhang, X., Hu, B., Qin, Y., Xu, H., Wang, H., Guo, X., Qian, J., et al. (2018).
# Root microbiota shift in rice correlates with resident time in the field and developmental stage. Sci China Life Sci 61, 
# https://doi.org/10.1007/s11427-018-9284-4

# 手动运行脚本请，需要设置工作目录，使用 Ctrl+Shift+H 或 Session - Set Work Directory - Choose Directory / To Source File Location 设置工作目录



# 1.1 程序功能描述和主要步骤

# 程序功能：高通量测序reads counts值的组间比较并筛选
# Functions: Calculate pvalue and FDR for each OTUs by edgeR or wilcon
# Main steps: 
# - Reads data matrix and design
# - Calculate pvalue and FDR
# - Save result table in *_all/sig.txt

# 清空工作环境 Clean enviroment object
rm(list=ls()) 


# 2.1 安装CRAN来源常用包
site="https://mirrors.tuna.tsinghua.edu.cn/CRAN"
# 依赖包列表：参数解析、数据变换、绘图和开发包安装、安装依赖、ggplot主题
package_list = c("limma","ggplot2","pheatmap","dplyr","devtools")
# 判断R包加载是否成功来决定是否安装后再加载
for(p in package_list){
	if(!suppressWarnings(suppressMessages(require(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
		install.packages(p, repos=site)
		suppressWarnings(suppressMessages(library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
  }
}

# 2.2 安装bioconductor常用包
package_list = c("edgeR")
for(p in package_list){
	if(!suppressWarnings(suppressMessages(require(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
		source("https://bioconductor.org/biocLite.R")
		biocLite(p)
		suppressWarnings(suppressMessages(library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
	}
}

# 2.3 安装Github常用包
# 参数解析、数据变换、绘图和开发包安装
package_list = c("kassambara/ggpubr")
for(p in package_list){
	q=unlist(strsplit(p,split = "/"))[2]
	if(!suppressWarnings(suppressMessages(require(q, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
		install_github(p)
		suppressWarnings(suppressMessages(library(q, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
	}
}


# 3. 读取输入文件

# 读取比较列表
input = read.table("result/compare/BacAhMH63wet-BacAhWYJ7DEP1wet_all.txt", header=T, row.names=1, sep="\t", comment.char="")
input$level=factor(input$level,levels = c("Enriched","Depleted","NotSig"))
NoE= dim(input[input$level=="Enriched",])[1]
NoD= dim(input[input$level=="Depleted",])[1]

# 绘制火山图，有上/下缺失时颜色会串
p = ggplot(input, aes(x=logFC, y=logCPM, color=level)) + 
	geom_point() + xlim(-4, 4) + theme_classic()+
	scale_colour_manual(values=c("red","green","grey")) + 
	labs(x="log2(fold change)", y="log2(count per million)")+ 
	annotate("text",x=-3,y=15,label=paste(NoD,sep=""))+ 
	annotate("text",x=3,y=15,label=paste(NoE,sep=""))
#p
suppressWarnings(ggsave(paste("result/compare/BacAhMH63wet-BacAhWYJ7DEP1wet", "_volcano.pdf", sep=""), p, width = 4, height = 2.5))
suppressWarnings(ggsave(paste("result/compare/BacAhMH63wet-BacAhWYJ7DEP1wet", "_volcano.png", sep=""), p, width = 4, height = 2.5))
# 提示工作完成
print(paste("Output in result/compare/BacAhMH63wet-BacAhWYJ7DEP1wet", "_volcano.pdf finished.", sep = ""))

