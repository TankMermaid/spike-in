# 2.1 安装CRAN来源常用包
site="https://mirrors.tuna.tsinghua.edu.cn/CRAN"
# 依赖包列表：参数解析、数据变换、绘图和开发包安装、安装依赖、ggplot主题
package_list = c("Biobase","edgeR","ggplot2","gplots","grid","RColorBrewer","reshape2","VennDiagram","dplyr","pheatmap")
# 判断R包加载是否成功来决定是否安装后再加载
for(p in package_list){
	if(!suppressWarnings(suppressMessages(require(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
		install.packages(p, repos=site)
		suppressWarnings(suppressMessages(library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
  }
}


# Set ggplot2 drawing parameter, such as axis line and text size, lengend and title size, and so on.
main_theme = theme(panel.background=element_blank(), panel.grid=element_blank(),
	axis.line.x=element_line(size=.5, colour="black"), axis.line.y=element_line(size=.5, colour="black"),
	axis.ticks=element_line(color="black"), axis.text=element_text(color="black", size=6),
	legend.position="right", legend.background=element_blank(), legend.key=element_blank(), legend.text= element_text(size=6),
	text=element_text(family="sans", size=6))

# 实验差异比较结果
x = read.table("result/compare/BacAhMH63wet-BacAhWYJ7DEP1wet_all.txt", header=T, row.names= 1, sep="\t", stringsAsFactors = F) 
# 只提取前14列
x = x[,1:14]
x = na.omit(x)

# P值求负对数
x$neglogp = -log10(x$PValue)

x$otu=rownames(x)
x = arrange(x, Kindom, Phylum, Class, Order, Family, Genus, Species)
x$otu = factor(x$otu, levels=x$otu)   # set x order
x$num = 1:dim(x)[1]

# 读取高丰度门，用于着色
# per= read.delim("result/tax/sum_p.txt", sep = "\t", row.names=1, header=T)
per= read.delim("result/tax/sum_pc.txt", sep = "\t", row.names=1, header=T)
mean = rowMeans(per)
per = as.data.frame(mean[order(mean, decreasing = T)])
top_tax=head(rownames(per), n=13)

# 将低丰度的门变为Low Abundance
# x$tax = factor(x$Phylum, levels=c(as.vector(unique(x$Phylum)),"Low Abundance"))

# 将低丰度的门变为Low Abundance
x$tax = x$Phylum# factor(x$Phylum, levels=c(as.vector(unique(x$Phylum)),"Low Abundance"))

# 将门中 proteobacteria 替换为纲
x[x$tax %in% "Proteobacteria",]$tax =  x[x$tax %in% "Proteobacteria",]$Class # no level can get value

if (length(unique(x$tax)) > length(top_tax)){
	x[!(x$tax %in% top_tax),]$tax = "Low Abundance" # no level can get value
}


# 调置标签顺序
label = unique(x$tax)
label = label[!(label %in% "Low Abundance")] # 删除低丰度
x$tax = factor(x$tax, levels = c(label, "Low Abundance"))
# 计算标签中位数
temp = x[x$tax %in% label, c("tax","num")]
mat_mean = aggregate(temp[,-1], by=temp[1], FUN=median) # mean


# 调整过大的负对数
if (max(x$neglogp)>20){
  x[x$neglogp>20,]$neglogp  = 20
}

# Manhattan plot
FDR = min(x$neglogp[x$level!="NotSig"])
p = ggplot(x, aes(x=num, y=neglogp, color=tax, size=logCPM, shape=level)) +
  geom_point(alpha=.7) + 
  geom_hline(yintercept=FDR, linetype=2, color="lightgrey") +
  scale_shape_manual(values=c(25, 17, 20))+
  scale_size(breaks=c(5, 10, 15)) +
  labs(x="Taxonomy of OTUs", y="-log10(P)", title=paste("result/compare/BacAhMH63wet-BacAhWYJ7DEP1wet_all.txt", sep=" ")) +main_theme +
  theme(legend.position="top") +
  scale_x_continuous(breaks=mat_mean$x, labels=mat_mean$tax) + 
  theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1))
p
ggsave(file=paste("result/compare/BacAhMH63wet-BacAhWYJ7DEP1wet_all.txt_man_pc.pdf", sep=""), p, width = 10, height = 4, useDingbats=F)
ggsave(file=paste("result/compare/BacAhMH63wet-BacAhWYJ7DEP1wet_all.txt_man_pc.png", sep=""), p, width = 10, height = 4)

