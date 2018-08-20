setwd("/mnt/bai/qinyuan/xiaoxuan/QA/")
print(paste("Your working directory is in",getwd()))
library(ggplot2)
library(tidyr)

## 读取数据集
data1 <- read.csv("data/Fig2b-bac-nospike-group2.csv", header = T)
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
ggsave(paste("result/Fig2b-bac-nospike-group2.pdf", sep=""), p, width = 8, height = 4)
