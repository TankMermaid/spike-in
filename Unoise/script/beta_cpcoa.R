#!/usr/bin/env Rscript
# 
# Copyright 2016-2018 Yong-Xin Liu <metagenome@126.com>

# If used this script, please cited:
#   Zhang, J., Zhang, N., Liu, Y.X., Zhang, X., Hu, B., Qin, Y., Xu, H., Wang, H., Guo, X., Qian, J., et al. (2018).
# Root microbiota shift in rice correlates with resident time in the field and developmental stage. Sci China Life Sci 61, 
# https://doi.org/10.1007/s11427-018-9284-4

# 手动运行脚本请，需要设置工作目录，使用 Ctrl+Shift+H 或 Session - Set Work Directory - Choose Directory / To Source File Location 设置工作目录



# 1.1 程序功能描述和主要步骤

# 程序功能：限制性主坐标轴分析及组间统计
# Functions: Constrained PCoA analysis of groups
# Main steps: 
# - Reads OTU table result/otutab.txt
# - Orrdinate by CCA and show in scatter plot
# - Aov.cca() calculate significant between all groups distance

# 清空工作环境 Clean enviroment object
rm(list=ls()) 


# 2.1 安装CRAN来源常用包
site="https://mirrors.tuna.tsinghua.edu.cn/CRAN"
# 依赖包列表：参数解析、数据变换、绘图和开发包安装、安装依赖、ggplot主题
package_list = c("reshape2","ggplot2","vegan")
# 判断R包加载是否成功来决定是否安装后再加载
for(p in package_list){
	if(!suppressWarnings(suppressMessages(require(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
		install.packages(p, repos=site)
		suppressWarnings(suppressMessages(library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
  }
}

# 2.2 安装bioconductor常用包
package_list = c("digest","ggrepel")
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

# Set ggplot2 drawing parameter, such as axis line and text size, lengend and title size, and so on.
main_theme = theme(panel.background=element_blank(), panel.grid=element_blank(),
	axis.line.x=element_line(size=.5, colour="black"), axis.line.y=element_line(size=.5, colour="black"),
	axis.ticks=element_line(color="black"), axis.text=element_text(color="black", size=7),
	legend.position="right", legend.background=element_blank(), legend.key=element_blank(), legend.text= element_text(size=7),
	text=element_text(family="sans", size=7))


# 3. 读取输入文件

# 读取实验设计
design = read.table("/mnt/bai/xiaoning/xiaoxuan/180528/180627_AQ/bac/unoise/doc/design.txt", header=T, row.names=1, sep="\t")
# 统一改实验列为group
design$group=design$groupID

# 按实验组筛选 Select by manual set group
if (TRUE){
	sub_design = subset(design, group %in% c("BacAhBulksoildry","BacAhBulksoilwet","BacAhMH63dry","BacAhMH63wet","BacAhMH63ZHdry","BacAhMH63ZHwet","BacAhWYJ7DEP1dry","BacAhWYJ7DEP1wet"))
# 调置组排序 Set group order
	sub_design$group  = factor(sub_design$group, levels=c("BacAhBulksoildry","BacAhBulksoilwet","BacAhMH63dry","BacAhMH63wet","BacAhMH63ZHdry","BacAhMH63ZHwet","BacAhWYJ7DEP1dry","BacAhWYJ7DEP1wet"))
}else{
	sub_design = design
}

# 函数1. 提取CCA中主要结果
# Function1. get CCA main result
variability_table = function(cca){
	chi = c(cca$tot.chi, cca$CCA$tot.chi, cca$CA$tot.chi)
	variability_table = cbind(chi, chi/chi[1])
	colnames(variability_table) = c("inertia", "proportion")
	rownames(variability_table) = c("total", "constrained", "unconstrained")
	return(variability_table)
}

# 4. 循环对每种距离矩阵统计和绘图

method = c("bray","jaccard")
for(m in method){
	# 读取usearch beta文件
	beta = read.table(paste("/mnt/bai/xiaoning/xiaoxuan/180528/180627_AQ/bac/unoise/result/otutab.txt",sep=""), header=T, row.names=1, sep="\t", comment.char="") 

	# 实验设计与输入文件交叉筛选
	idx = rownames(sub_design) %in% colnames(beta)
	sub_design=sub_design[idx,]
	sub_beta=beta[,rownames(sub_design)]

	# otutab标准化为10000
	otutab = t(sub_beta)/colSums(sub_beta,na=T)*10000

	# Constrained analysis OTU table by genotype
	capscale.gen = capscale(otutab ~ group, data=sub_design, add=F, sqrt.dist=T, distance= m) 

	# ANOVA-like permutation analysis
	perm_anova.gen = anova.cca(capscale.gen, permutations = 10000, parallel = 9)

	# generate variability tables and calculate confidence intervals for the variance
	var_tbl.gen = variability_table(capscale.gen)
	eig = capscale.gen$CCA$eig
	variance = var_tbl.gen["constrained", "proportion"]
	p.val = perm_anova.gen[1, 4]

	# extract the weighted average (sample) scores
	points = capscale.gen$CCA$wa[, 1:2]
	points = as.data.frame(points)
	points = cbind(points, sub_design$group)
	colnames(points) = c("PC1", "PC2", "group")

	# plot PC 1 and 2
	p = ggplot(points, aes(x=PC1, y=PC2, color=group)) + geom_point(alpha=.7, size=2) +
		labs(x=paste("CPCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
		y=paste("CPCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep="")) + 
		ggtitle(paste(format(100 * variance, digits=3), " % of variance; p=",format(p.val, digits=2),sep="")) + 
		theme_classic() + main_theme

	# 是否添加置信椭圆
	if (TRUE){
		p = p + stat_ellipse(level=0.68)
	}
	p
	# 保存pdf和png格式方便查看和编辑
	ggsave(paste("/mnt/bai/xiaoning/xiaoxuan/180528/180627_AQ/bac/unoise/result/beta/cpcoa_", m, ".pdf", sep=""), p, width = 8, height = 5)
	ggsave(paste("/mnt/bai/xiaoning/xiaoxuan/180528/180627_AQ/bac/unoise/result/beta/cpcoa_", m, ".png", sep=""), p, width = 8, height = 5)
	# 提示工作完成
	print(paste("Output in /mnt/bai/xiaoning/xiaoxuan/180528/180627_AQ/bac/unoise/result/beta/cpcoa_", m, ".pdf finished.", sep = ""))

	# 添加样品标签
	p = p + geom_text_repel(label = paste(rownames(points)), colour="black", size=3)
	p
	# 保存pdf和png格式方便查看和编辑
	# 保存pdf和png格式方便查看和编辑
	ggsave(paste("/mnt/bai/xiaoning/xiaoxuan/180528/180627_AQ/bac/unoise/result/beta/cpcoa_", m, "_label.pdf", sep=""), p, width = 8, height = 5)
#	ggsave(paste("/mnt/bai/xiaoning/xiaoxuan/180528/180627_AQ/bac/unoise/result/beta/cpcoa_", m, "_label.png", sep=""), p, width = 8, height = 5)
	# 提示工作完成
	print(paste("Output in /mnt/bai/xiaoning/xiaoxuan/180528/180627_AQ/bac/unoise/result/beta/cpcoa_", m, "_label.pdf finished.", sep = ""))
	print("")
}

