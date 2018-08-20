## 设置工作目录 加载包
setwd("/mnt/bai/qinyuan/xiaoxuan/bac_New/")
print(paste("Your working directory is in",getwd()))
library(ggplot2)
library(tidyr)

## 读取数据集
bacotu <- read.csv("table/Bac-ra.csv")
table <- gather(bacotu, key = spike , value = relativeabundance, `E00`:`E05.20`)

p <- ggplot(table, aes(spike, relativeabundance)) +
  geom_jitter(aes(colour = Bacteria), width = 0.08, na.rm = TRUE)+ 
  theme_classic()+
  scale_x_discrete(limits=c("E00","E04mix","E05","E05.5","E05.10","E05.20"))
p
ggsave(paste("figure/bacwithoutspike.pdf", sep=""), p, width = 6, height = 5)
