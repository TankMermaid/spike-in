#!/usr/bin/env Rscript
# 
# Copyright 2016-2018 Yong-Xin Liu <metagenome@126.com>

# If used this script, please cited:
#   Zhang, J., Zhang, N., Liu, Y.X., Zhang, X., Hu, B., Qin, Y., Xu, H., Wang, H., Guo, X., Qian, J., et al. (2018).
# Root microbiota shift in rice correlates with resident time in the field and developmental stage. Sci China Life Sci 61, 
# https://doi.org/10.1007/s11427-018-9284-4

# 手动运行脚本请，需要设置工作目录，使用 Ctrl+Shift+H 或 Session - Set Work Directory - Choose Directory / To Source File Location 设置工作目录



# 1.1 程序功能描述和主要步骤

# 程序功能：Beta多样性主坐标轴分析及组间统计
# Functions: PCoA analysis of samples and groups comparing
# Main steps: 
# - Reads distance matrix input.txt
# - Calculate orrdinate by PCoA and show in scatter plot
# - Adonis calculate significant between groups distance and group inner distance

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

# 函数1：采用adonis对距离矩阵进行组间差异统计
# Compare each group distance matrix by vegan adonis in bray_curtis
da_adonis = function(sampleV){
	sampleA = as.matrix(sampleV$sampA)
	sampleB = as.matrix(sampleV$sampB)
	design2 = subset(sub_design, group %in% c(sampleA,sampleB))
	if (length(unique(design2$group))>1) {
		sub_dis_table = dis_table[rownames(design2),rownames(design2)]
		sub_dis_table = as.dist(sub_dis_table, diag = FALSE, upper = FALSE)
		adonis_table = adonis(sub_dis_table ~ group, data = design2, permutations = 1000) 
		adonis_pvalue = adonis_table$aov.tab$`Pr(>F)`[1]
		print(paste("In", m, "pvalue between", sampleA, "and", sampleB, "is", adonis_pvalue, sep=" "))
		adonis_pvalue = paste(m, sampleA, sampleB, adonis_pvalue, sep="\t")
		return(adonis_pvalue)
	}
}


# 3. 读取输入文件

# 读取实验设计
design = read.table("/mnt/bai/yongxin/other/guoxiaoxuan/180627_AQ/bac/doc/design.txt", header=T, row.names=1, sep="\t")
# 统一改实验列为group
design$group=design$groupID

# 按实验组筛选 Select by manual set group
if (TRUE){
	sub_design = subset(design, group %in% c("BacHnBulksoildry","BacHnBulksoilwet","BacHnMH63dry","BacHnMH63wet","BacHnMH63ZHdry","BacHnMH63ZHwet","BacHnWYJ7DEP1dry","BacHnWYJ7DEP1wet"))
# 调置组排序 Set group order
	sub_design$group  = factor(sub_design$group, levels=c("BacHnBulksoildry","BacHnBulksoilwet","BacHnMH63dry","BacHnMH63wet","BacHnMH63ZHdry","BacHnMH63ZHwet","BacHnWYJ7DEP1dry","BacHnWYJ7DEP1wet"))
}else{
	sub_design = design
}



# 4. 循环对每种距离矩阵统计和绘图


method = c("bray_curtis","weighted_unifrac","unweighted_unifrac")
for(m in method){
	# 读取usearch beta文件
	beta = read.table(paste("/mnt/bai/yongxin/other/guoxiaoxuan/180627_AQ/bac/spikein/beta/",m,".txt",sep=""), header=T, row.names=1, sep="\t", comment.char="") 
	write.table(paste("#Permutational Multivariate Analysis of Variance Using Distance Matrices\nDistanceMatrices\tGroupA\tGroupB\tP-value",  sep=""), 
				file=paste("/mnt/bai/yongxin/other/guoxiaoxuan/180627_AQ/bac/spikein/beta/Hn", m, ".stat", sep=""), append = F, sep="\t", quote=F, row.names=F, col.names=F)

	# 实验设计与输入文件交叉筛选
	idx = rownames(sub_design) %in% rownames(beta)
	sub_design=sub_design[idx,]
	sub_beta=beta[rownames(sub_design),rownames(sub_design)]

	# vegan:cmdscale计算矩阵矩阵中主坐标轴坐标，取前3维
	pcoa = cmdscale(sub_beta, k=4, eig=T) # k is dimension, 3 is recommended; eig is eigenvalues
	points = as.data.frame(pcoa$points) # get coordinate string, format to dataframme
	eig = pcoa$eig
	points = cbind(points, sub_design$group)
	colnames(points) = c("PC1", "PC2", "PC3", "PC4","group") 
	write.table("Samples\t", file=paste("/mnt/bai/yongxin/other/guoxiaoxuan/180627_AQ/bac/spikein/beta/Hn", m, "14.txt",sep = ""), append = F, sep="\t", quote=F,  eol = "",row.names=F, col.names=F)
	suppressWarnings(write.table(points[,c(1:4)], file=paste("/mnt/bai/yongxin/other/guoxiaoxuan/180627_AQ/bac/spikein/beta/Hn", m, "14.txt",sep = ""), append = T, sep="\t", quote=F, row.names=T, col.names=T))

	# plot PC 1 and 2
	p = ggplot(points, aes(x=PC1, y=PC2, color=group)) + geom_point(alpha=.7, size=2) +
		labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
		y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep=""),
		title=paste(m," PCoA",sep="")) + theme_classic()

	# 是否添加置信椭圆
	if (TRUE){
		p = p + stat_ellipse(level=0.68)
	}
	p
	# 保存pdf和png格式方便查看和编辑
	ggsave(paste("/mnt/bai/yongxin/other/guoxiaoxuan/180627_AQ/bac/spikein/beta/Hn", m, ".pdf", sep=""), p, width = 8, height = 5)
	ggsave(paste("/mnt/bai/yongxin/other/guoxiaoxuan/180627_AQ/bac/spikein/beta/Hn", m, ".png", sep=""), p, width = 8, height = 5)
	# 提示工作完成
	print(paste("Output in /mnt/bai/yongxin/other/guoxiaoxuan/180627_AQ/bac/spikein/beta/Hn", m, ".pdf finished.", sep = ""))

	# 添加样品标签
	p = p + geom_text_repel(label = paste(rownames(points)), colour="black", size=3)
	p
	# 保存pdf和png格式方便查看和编辑
	# 保存pdf和png格式方便查看和编辑
	ggsave(paste("/mnt/bai/yongxin/other/guoxiaoxuan/180627_AQ/bac/spikein/beta/Hn", m, "_label.pdf", sep=""), p, width = 8, height = 5)
#	ggsave(paste("/mnt/bai/yongxin/other/guoxiaoxuan/180627_AQ/bac/spikein/beta/Hn", m, "_label.png", sep=""), p, width = 8, height = 5)
	# 提示工作完成
	print(paste("Output in /mnt/bai/yongxin/other/guoxiaoxuan/180627_AQ/bac/spikein/beta/Hn", m, "_label.pdf finished.", sep = ""))

	# plot PC 3 and 4
	p = ggplot(points, aes(x=PC3, y=PC4, color=group)) + geom_point(alpha=.7, size=2) +
		labs(x=paste("PCoA 3 (", format(100 * eig[3] / sum(eig), digits=4), "%)", sep=""),
		y=paste("PCoA 4 (", format(100 * eig[4] / sum(eig), digits=4), "%)", sep=""),
		title=paste(m," PCoA",sep="")) + theme_classic()
	# 是否添加置信椭圆
	if (TRUE){p = p + stat_ellipse(level=0.68)}
	p
	# 保存pdf和png格式方便查看和编辑
	ggsave(paste("/mnt/bai/yongxin/other/guoxiaoxuan/180627_AQ/bac/spikein/beta/Hn", m, "_34.pdf", sep=""), p, width = 8, height = 5)
#	ggsave(paste("/mnt/bai/yongxin/other/guoxiaoxuan/180627_AQ/bac/spikein/beta/Hn", m, "_34.png", sep=""), p, width = 8, height = 5)
	# 提示工作完成
	print(paste("Output in /mnt/bai/yongxin/other/guoxiaoxuan/180627_AQ/bac/spikein/beta/Hn", m, "_34.pdf finished.", sep = ""))

	# 循环统计各比较组或所有组 loop for each group pair
	dis_table = as.matrix(sub_beta)
	# 如果没有比较文件，则自动全循环
	if (!file.exists("/mnt/bai/yongxin/other/guoxiaoxuan/180627_AQ/bac/doc/compare.txt")) {
		compare_data = as.vector(unique(sub_design$group))
		len_compare_data = length(compare_data)
		for(i in 1:(len_compare_data-1)) {
			for(j in (i+1):len_compare_data) {
				tmp_compare = as.data.frame(cbind(sampA=compare_data[i],sampB=compare_data[j]))
#				print(tmp_compare)
				adonis_pvalue = da_adonis(tmp_compare)
				write.table(adonis_pvalue, file=paste("/mnt/bai/yongxin/other/guoxiaoxuan/180627_AQ/bac/spikein/beta/Hn", m, ".stat", sep=""), append = TRUE, sep="\t", quote=F, row.names=F, col.names=F)
			}
		}
	# 有比较文件，按设计比较
	}else {
		compare_data = read.table("/mnt/bai/yongxin/other/guoxiaoxuan/180627_AQ/bac/doc/compare.txt", sep="\t", check.names=F, quote='', comment.char="")
		colnames(compare_data) = c("sampA", "sampB")
		for(i in 1:dim(compare_data)[1]){
			adonis_pvalue = da_adonis(compare_data[i,])
			write.table(adonis_pvalue, file=paste("/mnt/bai/yongxin/other/guoxiaoxuan/180627_AQ/bac/spikein/beta/Hn", m, ".stat", sep=""), append = TRUE, sep="\t", quote=F, row.names=F, col.names=F)
		}
	}	 
	print(paste("Adnois statistics result in ", "/mnt/bai/yongxin/other/guoxiaoxuan/180627_AQ/bac/spikein/beta/Hn", m, ".stat\n is finished.", sep = ""))
	print("")
}

