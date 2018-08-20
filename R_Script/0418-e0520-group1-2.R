## 设置工作目录 加载包
setwd("/mnt/bai/qinyuan/xiaoxuan/bac_New/")
print(paste("Your working directory is in",getwd()))
library(ggplot2)
library(tidyr)

## 读取数据集
rawtable <- read.csv("table/20180418-e0520-group1-2.csv", header = T)
colnames(rawtable)[1] = "bacteria"
table <- gather(rawtable, key = sample , value = count, `sample01`:`sample05`)

## 设置横轴变量为离散型因子
table1 <- table
table1$bacteria <- factor(table1$bacteria)
table1$group <- factor(table1$group)

## 生成箱线图
boxplot1 <- ggplot(table1) +
  geom_boxplot(aes(x = bacteria, y = count, fill = group), width = 0.3, na.rm = TRUE) +
  geom_line(aes(x = bacteria, y = mean, group = group, color = group), na.rm = TRUE) +
  labs(x="Bacteria", y="Relative abundance (%)") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust = 0))
boxplot1
ggsave(paste("figure/20180418-e0520-group1-2.pdf", sep=""), boxplot1, width = 6, height = 5)
