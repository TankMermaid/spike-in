setwd("/mnt/bai/qinyuan/xiaoxuan/QA/")
print(paste("Your working directory is in",getwd()))
rm(list=ls()) 

library(ggplot2)
library(ggpubr)
library(dplyr)

data = read.table("./data/Fig6e.txt",header = T, row.names = 1)
head(data)

data_ra <- filter(data, method == "ra")

p <- ggboxplot(data_ra, x="treatment",
               y = "count", color = "treatment",
               palette =  "npg", add = "jitter",
               facet.by = "genus", short.panel.labs = FALSE)
# 添加p值
p01 <- p + stat_compare_means(aes(label = ..p.signif..), method = "t.test")

p01
ggsave(paste("./result/Fig6e-1.pdf", sep=""), p01, width = 5, height = 8)


data_aa <- filter(data, method == "aa")

p2 <- ggboxplot(data_aa, x="treatment",
               y = "count", color = "treatment",
               palette = "jco", add = "jitter",
               facet.by = "genus", short.panel.labs = FALSE)
# 添加p值
p02 <- p2 + stat_compare_means(aes(label = ..p.signif..),method = "t.test")

p02
ggsave(paste("./result/Fig6e-2.pdf", sep=""), p02, width = 5, height = 8)