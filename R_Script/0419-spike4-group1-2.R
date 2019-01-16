## 设置工作目录 加载包
setwd("/mnt/bai/qinyuan/xiaoxuan/bac_New/")
print(paste("Your working directory is in",getwd()))
library(ggplot2)
library(tidyr)

## 读取数据集
rawtable <- read.csv("table/20180419-spike4-group1-2.csv", header = T)
colnames(rawtable)[1] = "bacteria"
table <- gather(rawtable, key = sample , value = count, `sample01`:`sample05`)

## 设置横轴变量为离散型因子
table3 <- table
table3$bacteria <- factor(table$bacteria)
table3$group <- factor(table$group)

## 生成箱线图
boxplot3 <- ggplot(table3) +
  geom_boxplot(aes(x = bacteria, y = count, fill = group), width = 0.3, na.rm = TRUE) +
  geom_line(aes(x = bacteria, y = mean, group = group, color = group), na.rm = TRUE) +
  facet_wrap(~spike, ncol = 4) +
  labs(x="Bacteria", y="Relative abundance (%)") +
  theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust = 0))
boxplot3

ggsave(paste("figure/20180419-spike4-group1-2.pdf", sep=""), boxplot3, width = 9, height = 5)
