setwd("/mnt/bai/qinyuan/xiaoxuan/bac_New/")
print(paste("Your working directory is in",getwd()))
library(ggplot2)
library(tidyr)

## 读取数据集
data1 <- read.csv("table/20180413-bac-nospike-01-05.csv", header = T)
colnames(data1)[1] = "bacteria"

data2 <- gather(data1, key = replicate , value = count, `sample01`:`sample05`)

## 作图
p <- ggplot(data2, aes(bacteria, count))+
  geom_boxplot(aes(fill = bacteria), na.rm = TRUE) +
  facet_wrap(~spike, ncol = 5) +
  scale_fill_manual(values = c("#ff8364", "#fdb87d", "#ffe8d5", "#fbd685","#7fa99b", 
                               "#d5eeff", "#87e0ff", "#53c7f0", "#1d97c1"))+
  labs(x="Bacteria", y="Count (%)")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust = 0)) 
p
ggsave(paste("figure/20180413-bac-nospike-01-05.pdf", sep=""), p, width = 8, height = 4)

## 读取数据集2
data3 <- read.csv("table/20180413-bac-nospike-06-10.csv", header = T)
colnames(data3)[1] = "bacteria"

data4 <- gather(data3, key = replicate , value = count, `sample06`:`sample10`)

## 作图
p2 <- ggplot(data4, aes(bacteria, count)) + 
  geom_boxplot(aes(fill = bacteria), na.rm = TRUE) +
  facet_wrap(~spike, ncol = 5)  +
  scale_fill_manual(values = c("#ff8364", "#fdb87d", "#ffe8d5", "#fbd685","#7fa99b", 
                               "#d5eeff", "#87e0ff", "#53c7f0", "#1d97c1"))+
  labs(x="Bacteria", y="Count (%)")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust = 0)) 
p2
ggsave(paste("figure/20180413-bac-nospike-06-10.pdf", sep=""), p2, width = 8, height = 4)

## 读取数据集
data5 <- read.csv("table/20180413-bac-nospike-11-15.csv", header = T)
colnames(data5)[1] = "bacteria"

data6 <- gather(data5, key = replicate , value = count, `sample11`:`sample15`)

## 作图
p3 <- ggplot(data6, aes(bacteria, count)) + 
  geom_boxplot(aes(fill = bacteria), na.rm = TRUE) +
  facet_wrap(~spike, ncol = 5) +
  theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust = 0)) +
  scale_fill_brewer(palette = "Set3") +
  labs(x="Bacteria", y="Count (%)")
p3
ggsave(paste("figure/20180413-bac-nospike-11-15.pdf", sep=""), p3, width = 8, height = 4)


