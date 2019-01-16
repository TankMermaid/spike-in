## 设置工作目录 加载包
setwd("/mnt/bai/qinyuan/xiaoxuan/bac_New/")
print(paste("Your working directory is in",getwd()))
library(ggplot2)
library(tidyr)

## 读取数据集
group <- read.csv("table/20180417-ra-e0520.csv", header = T)
colnames(group)[1] = "bacteria"

## 设置横轴变量为离散型因子
dataset <- group
dataset$bacteria <- factor(dataset$bacteria)


## 作图

p <- ggplot(dataset) +
  geom_line(aes(x = bacteria, y = mean1, group = 1, color = "2:2:2")) +
  geom_line(aes(x = bacteria, y = mean2, group = 1, color = "1:1:1"))
p


p2 <- ggplot(dataset) +
  geom_line(aes(x = bacteria, y = mean1, group = 1, colour = "2:2:2")) +
  geom_line(aes(x = bacteria, y = mean3, group = 1, colour = "2:2:1"))
p2
