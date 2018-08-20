setwd("/mnt/bai/qinyuan/xiaoxuan/QA/")
print(paste("Your working directory is in",getwd()))
library(ggplot2)
library(dplyr)

load <- read.csv("data/Fig5f.csv", header = T)

## AH 竖版

data_AH <- filter(load, Site == "AH")

p3 <- ggplot(data_AH , aes(x= Type, y = Mean, fill = Micro )) + 
  geom_bar(stat = "identity", width=0.7) + 
  geom_errorbar(aes(ymin=Mean2-Error, ymax=Mean2+Error), width=.1) +
  labs(x="Group", y="Quantitative kingdom abundance (QA)") +
  facet_wrap(~Genotype, ncol = 2) +
  scale_fill_manual(values = c("#fdffab","#a8e6cf"))+
  theme_bw()
p3
ggsave(paste("result/Fig5f-data_AH.pdf", sep=""), p3, width = 5, height = 6)


## HN

data_HN <- filter(load, Site == "HN")

p4 <- ggplot(data_HN , aes(x= Type, y = Mean, fill = Micro )) + 
  geom_bar(stat = "identity", width=0.7) + 
  geom_errorbar(aes(ymin=Mean2-Error, ymax=Mean2+Error), width=.1) +
  labs(x="Group", y="Quantitative kingdom abundance (QA)") +
  facet_wrap(~Genotype, ncol = 2) +
  scale_fill_manual(values = c("#fdffab","#a8e6cf"))+
  theme_bw()
p4
ggsave(paste("result/Fig5f-data_HN.pdf", sep=""), p4, width = 5, height = 6)


