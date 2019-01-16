## 设置工作目录 加载包
setwd("/mnt/bai/qinyuan/xiaoxuan/bac_New/")
print(paste("Your working directory is in",getwd()))
library(ggplot2)
library(tidyr)

## 读取数据集
rawtable1 <- read.csv("table/20180419-spike4-group1-2.csv", header = T)
colnames(rawtable1)[1] = "bacteria"
ra1 <- gather(rawtable1, key = sample , value = count, `sample01`:`sample05`)

rawtable2 <- read.csv("table/20180419-spike4-group1-3.csv", header = T)
colnames(rawtable2)[1] = "bacteria"
ra2 <- gather(rawtable2, key = sample , value = count, `sample01`:`sample05`)

rawtable3 <- read.csv("table/20180420-spike3-group1-2.csv", header = T)
colnames(rawtable3)[1] = "bacteria"
aa1 <- gather(rawtable3, key = sample , value = count, `sample01`:`sample05`)

rawtable4 <- read.csv("table/20180420-spike3-group1-3.csv", header = T)
colnames(rawtable4)[1] = "bacteria"
aa2 <- gather(rawtable4, key = sample , value = count, `sample01`:`sample05`)

## 设置横轴变量为离散型因子
ra1data <- ra1
ra1data$bacteria <- factor(ra1data$bacteria)
ra1data$group <- factor(ra1data$group)

ra2data <- ra2
ra2data$bacteria <- factor(ra2data$bacteria)
ra2data$group <- factor(ra2data$group)

aa1data <- aa1
aa1data$bacteria <- factor(aa1data$bacteria)
aa1data$group <- factor(aa1data$group)

aa2data <- aa2
aa2data$bacteria <- factor(aa2data$bacteria)
aa2data$group <- factor(aa2data$group)

## 生成折线图
line1 <- ggplot(ra1data) +
  geom_point(aes(x = bacteria, y = mean, color = group), na.rm = TRUE) +
  geom_line(aes(x = bacteria, y = mean, group = group, color = group), na.rm = TRUE) +
  facet_wrap(~spike, ncol = 4) +
  labs(x="Bacteria", y="Relative abundance (%)") +
  theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust = 0))
line1

line2 <- ggplot(ra2data) +
  geom_point(aes(x = bacteria, y = mean, color = group), na.rm = TRUE) +
  geom_line(aes(x = bacteria, y = mean, group = group, color = group), na.rm = TRUE) +
  facet_wrap(~spike, ncol = 4) +
  labs(x="Bacteria", y="Relative abundance (%)") +
  theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust = 0))
line2

line3 <- ggplot(aa1data) +
  geom_point(aes(x = bacteria, y = mean, color = group), na.rm = TRUE) +
  geom_line(aes(x = bacteria, y = mean, group = group, color = group), na.rm = TRUE) +
  facet_wrap(~spike, ncol = 4) +
  labs(x="Bacteria", y="Absolute abundance (%)") +
  theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust = 0))
line3

line4 <- ggplot(aa2data) +
  geom_point(aes(x = bacteria, y = mean, color = group), na.rm = TRUE) +
  geom_line(aes(x = bacteria, y = mean, group = group, color = group), na.rm = TRUE) +
  facet_wrap(~spike, ncol = 4) +
  labs(x="Bacteria", y="Absolute abundance (%)") +
  theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust = 0))
line4

## 保存
ggsave(paste("figure/20180420-line-facet1.pdf", sep=""), line1, width = 10, height = 5)
ggsave(paste("figure/20180420-line-facet2.pdf", sep=""), line2, width = 10, height = 5)
ggsave(paste("figure/20180420-line-facet3.pdf", sep=""), line3, width = 10, height = 5)
ggsave(paste("figure/20180420-line-facet4.pdf", sep=""), line4, width = 10, height = 5)


