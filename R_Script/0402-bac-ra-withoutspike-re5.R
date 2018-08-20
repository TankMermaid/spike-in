## 设置工作目录 加载包
setwd("/mnt/bai/qinyuan/xiaoxuan/bac_New/")
print(paste("Your working directory is in",getwd()))
library(ggplot2)
library(tidyr)

## 读取和整理数据集
bacotu <- read.csv("table/GOS.Bac.re5.csv", header = T)
table <- data.frame(bacotu$Bacteria, bacotu$spike, bacotu$averge, bacotu$stdev,row.names = row.names(bacotu))
colnames(table)[1] = "bacteria"
colnames(table)[2] = "spike"
colnames(table)[3] = "average"
colnames(table)[4] = "stdev"

p <- ggplot(table, aes(spike, average, ymin = average - stdev, ymax = average + stdev, colour = bacteria)) +
  geom_pointrange(position = position_jitter(width = 0.08)) +
  theme_classic()+
  scale_x_discrete(limits=c("E00","E04mix","E05","E05-5","E05-10","E05-20"))
p
ggsave(paste("figure/bacwithoutspikere5.pdf", sep=""), p, width = 6, height = 5)
