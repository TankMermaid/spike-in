## 设置工作目录 加载包
setwd("/mnt/bai/qinyuan/xiaoxuan/bac_New/")
print(paste("Your working directory is in",getwd()))
library(ggplot2)
library(tidyr)

## 读取数据集
group <- read.csv("table/20180418-e0520.csv", header = T)
colnames(group)[1] = "bacteria"
data <- gather(group, key = sample , value = count, `sample01`:`sample05`)

## 设置横轴变量为离散型因子
dataset <- data
dataset$bacteria <- factor(dataset$bacteria)
dataset$group <- factor(dataset$group)


p <- ggplot(dataset) +
  geom_boxplot(aes(x = bacteria, y = count, fill = group), width = 0.3, na.rm = TRUE) +
  geom_line(aes(x = bacteria, y = mean, group = 1, colour = group)) +
  labs(x="Bacteria", y="Relative abundance (%)") +
  theme_classic()
p



