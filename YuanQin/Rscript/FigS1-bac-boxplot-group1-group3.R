setwd("/mnt/bai/qinyuan/xiaoxuan/QA/")
print(paste("Your working directory is in",getwd()))
library(ggplot2)
library(tidyr)
library(dplyr)

## 读取数据集
data <- read.csv("data/FigS1-bac-nospike-group1-group3.csv", header = T)
colnames(data)[1] = "bacteria"
group1 <- select(data, bacteria, spike, sample06, sample09, sample10)
group3 <- select(data, bacteria, spike, sample11, sample12, sample13, sample14, sample15)

## group 1

data1 <- gather(group1, key = replicate , value = count, `sample06`:`sample10`)

p1 <- ggplot(data1, aes(bacteria, count))+
  geom_boxplot(aes(fill = bacteria), na.rm = TRUE) +
  facet_wrap(~spike, ncol = 5) +
  scale_fill_manual(values = c("#ff8364", "#fdb87d", "#ffe8d5", "#fbd685","#7fa99b", 
                               "#d5eeff", "#87e0ff", "#53c7f0", "#1d97c1"))+
  labs(x="Bacteria", y="Count (%)")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust = 0)) 
p1

ggsave(paste("result/FigS1-bac-group1.pdf", sep=""), p1, width = 8, height = 4)


## group 3

data3 <- gather(group3, key = replicate , value = count, `sample11`:`sample15`)

p2 <- ggplot(data3, aes(bacteria, count))+
  geom_boxplot(aes(fill = bacteria), na.rm = TRUE) +
  facet_wrap(~spike, ncol = 5) +
  scale_fill_manual(values = c("#ff8364", "#fdb87d", "#ffe8d5", "#fbd685","#7fa99b", 
                               "#d5eeff", "#87e0ff", "#53c7f0", "#1d97c1"))+
  labs(x="Bacteria", y="Count (%)")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust = 0)) 
p2

ggsave(paste("result/FigS1-bac-group3.pdf", sep=""), p2, width = 8, height = 4)

