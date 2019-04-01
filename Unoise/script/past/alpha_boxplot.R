#!/usr/bin/env Rscript
# 
# Copyright 2016-2018 Yong-Xin Liu <metagenome@126.com>

# If used this script, please cited:
#   Zhang, J., Zhang, N., Liu, Y.X., Zhang, X., Hu, B., Qin, Y., Xu, H., Wang, H., Guo, X., Qian, J., et al. (2018).
# Root microbiota shift in rice correlates with resident time in the field and developmental stage. Sci China Life Sci 61, 
# https://doi.org/10.1007/s11427-018-9284-4

# 手动运行脚本请，需要设置工作目录，使用 Ctrl+Shift+H 或 Session - Set Work Directory - Choose Directory / To Source File Location 设置工作目录



# 1.1 程序功能描述和主要步骤

# 程序功能：Alpha多样性箱线图绘制及组间统计
# Functions: Boxplot show alpha richness among groups
# Main steps: 
# - Reads data table input.txt
# - Calculate pvalue and save in output.txt
# - Draw boxplot and save in output.pdf



# 2.1 安装CRAN来源常用包
site="https://mirrors.tuna.tsinghua.edu.cn/CRAN"
# 依赖包列表：参数解析、数据变换、绘图和开发包安装、安装依赖、ggplot主题
package_list = c("reshape2","ggplot2","devtools","bindrcpp",
				"ggthemes","agricolae","dplyr")
# 判断R包加载是否成功来决定是否安装后再加载
for(p in package_list){
	if(!suppressWarnings(suppressMessages(require(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
		install.packages(p, repos=site)
		suppressWarnings(suppressMessages(library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
  }
}

# 2.2 安装bioconductor常用包
package_list = c("digest")
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
	axis.ticks=element_line(color="black"), axis.text=element_text(color="black", size=),
	legend.position="right", legend.background=element_blank(), legend.key=element_blank(), legend.text= element_text(size=),
	text=element_text(family="sans", size=))


# 3. 读取输入文件

# 读取usearch alpha文件
alpha = read.table("/mnt/bai/yongxin/other/guoxiaoxuan/180627_AQ/bac/result/alpha/index.txt", header=T, row.names=1, sep="\t", comment.char="") 

# 转置数据矩阵
if (FALSE){
	alpha = as.data.frame(t(alpha))
}

# 标准化数据矩阵
if (FALSE){
	alpha = as.data.frame(alpha/rowSums(alpha,na=T) * 100) # normalization to 1000
}else{
	alpha = alpha * 1
}


# 读取实验设计
design = read.table("/mnt/bai/yongxin/other/guoxiaoxuan/180627_AQ/bac/doc/design.txt", header=T, row.names=1, sep="\t")
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

# 实验设计与输入文件交叉筛选
idx = rownames(sub_design) %in% rownames(alpha)
sub_design=sub_design[idx, , drop=F]
sub_alpha=alpha[rownames(sub_design),]

# 合并Alpha指数与实验设计 add design to alpha
index = cbind(sub_alpha, sub_design) 



# 4. 循环对每种指数统计和绘图
method = c("chao1","richness","shannon_e")
for(m in method){
	# 统计各组间差异
	model = aov(index[[m]] ~ group, data=index)
	# 计算Tukey显著性差异检验
	Tukey_HSD = TukeyHSD(model, ordered = TRUE, conf.level = 0.95)
	# 提取比较结果
	Tukey_HSD_table = as.data.frame(Tukey_HSD$group) 
	# 保存一个制表符，解决存在行名时，列名无法对齐的问题
	write.table(paste(m, "\n\t", sep=""), file=paste("/mnt/bai/yongxin/other/guoxiaoxuan/180627_AQ/bac/result/alpha/",m,".txt",sep=""),append = F, quote = F, eol = "", row.names = F, col.names = F)
	# 保存统计结果，有waring正常
	suppressWarnings(write.table(Tukey_HSD_table, file=paste("/mnt/bai/yongxin/other/guoxiaoxuan/180627_AQ/bac/result/alpha/",m,".txt",sep=""), append = T, quote = F, sep="\t", eol = "\n", na = "NA", dec = ".", row.names = T, col.names = T))

	# LSD检验，添加差异组字母
	out = LSD.test(model,"group", p.adj="none") # alternative fdr
	stat = out$groups
	# 分组结果添入Index
	index$stat=stat[as.character(index$group),]$groups
	# 设置分组位置为各组y最大值+高的5%
	max=max(index[,c(m)])
	min=min(index[,c(m)])
	x = index[,c("group",m)]
	y = x %>% group_by(group) %>% summarise_(Max=paste('max(',m,')',sep=""))
	y=as.data.frame(y)
	rownames(y)=y$group
	index$y=y[as.character(index$group),]$Max + (max-min)*0.05
	
# 输出原始数据，方便筛选 
if (FALSE){
	write.table(paste("SampleID\t", sep=""), file=paste("/mnt/bai/yongxin/other/guoxiaoxuan/180627_AQ/bac/result/alpha/",m,"_raw.txt",sep=""), append = F, quote = F, eol = "", row.names = F, col.names = F)
	suppressWarnings(write.table(index[,c(m,"group")], file=paste("/mnt/bai/yongxin/other/guoxiaoxuan/180627_AQ/bac/result/alpha/",m,"_raw.txt",sep=""), append = T, quote = F, sep="\t", eol = "\n", na = "NA", dec = ".", row.names = T, col.names = T))
}

	p = ggplot(index, aes(x=group, y=index[[m]], color=group)) +
		geom_boxplot(alpha=1, outlier.size=0, size=0.7, width=0.5, fill="transparent") +
		labs(x="Groups", y=paste(m, "")) + theme_classic() + main_theme +
		geom_text(data=index, aes(x=group, y=y, color=group, label= stat)) +
		geom_jitter( position=position_jitter(0.17), size=1, alpha=0.7)
	if (length(unique(sub_design$group))>3){
		p=p+theme(axis.text.x=element_text(angle=45,vjust=1, hjust=1))
	}


	p
	ggsave(paste("/mnt/bai/yongxin/other/guoxiaoxuan/180627_AQ/bac/result/alpha/", m, ".pdf", sep=""), p, width = 8, height = 5)
	ggsave(paste("/mnt/bai/yongxin/other/guoxiaoxuan/180627_AQ/bac/result/alpha/", m, ".png", sep=""), p, width = 8, height = 5)
	# 提示工作完成
	print(paste("Output in /mnt/bai/yongxin/other/guoxiaoxuan/180627_AQ/bac/result/alpha/", m, ".txt/pdf finished.", sep = ""))
}

