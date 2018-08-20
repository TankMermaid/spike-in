## 设置工作目录 加载包
setwd("/mnt/bai/qinyuan/xiaoxuan/bac_New/")
print(paste("Your working directory is in",getwd()))
library(ggplot2)
library(tidyr)

## 读取数据集
rawtable <- read.csv("table/20180420-e0520-group1-2.csv", header = T)
colnames(rawtable)[1] = "bacteria"
aatable <- gather(rawtable, key = sample , value = count, `sample01`:`sample05`)

## 设置横轴变量为离散型因子
aatable1 <- aatable
aatable1$bacteria <- factor(aatable1$bacteria)
aatable1$group <- factor(aatable1$group)

## 生成箱线图
aaboxplot1 <- ggplot(aatable1) +
  geom_boxplot(aes(x = bacteria, y = count, fill = group), width = 0.3, na.rm = TRUE) +
  geom_line(aes(x = bacteria, y = mean, group = group, color = group), na.rm = TRUE) +
  labs(x="Bacteria", y="Absolute abundance (%)") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust = 0))
aaboxplot1
ggsave(paste("figure/20180420-e0520-group1-2.pdf", sep=""), aaboxplot1, width = 6, height = 5)
