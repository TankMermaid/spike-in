setwd("/mnt/bai/qinyuan/xiaoxuan/bac_New/")
print(paste("Your working directory is in",getwd()))
library(ggplot2)
library(dplyr)

load <- read.csv("table/20180804-FIG5f.csv", header = T)

## AH 横版

data_AH <- filter(load, Site == "AH")

p1 <- ggplot(data_AH , aes(x= Type, y = Mean, fill = Micro )) + 
  geom_bar(stat = "identity", width=0.7) + 
  geom_errorbar(aes(ymin=Mean2-Error, ymax=Mean2+Error), width=.1) +
  labs(x="Group", y="Quantitative kingdom abundance (QA)") +
  facet_grid(Genotype ~ ., scales = "free", space = "free") +
  scale_fill_manual(values = c("#fdffab","#a8e6cf"))+
  coord_flip()+
  theme_bw()
p1
ggsave(paste("figure/20180804-fig5f-data_AH.pdf", sep=""), p1, width = 7, height = 4)


## HN

data_HN <- filter(load, Site == "HN")

p2 <- ggplot(data_HN , aes(x= Type, y = Mean, fill = Micro )) + 
  geom_bar(stat = "identity", width=0.7) + 
  geom_errorbar(aes(ymin=Mean2-Error, ymax=Mean2+Error), width=.1) +
  labs(x="Group", y="Quantitative kingdom abundance (QA)") +
  facet_grid(Genotype ~ ., scales = "free", space = "free") +
  scale_fill_manual(values = c("#fdffab","#a8e6cf"))+
  coord_flip()+
  theme_bw()
p2
ggsave(paste("figure/20180804-fig5f-data_HN.pdf", sep=""), p2, width = 7, height = 4)


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
ggsave(paste("figure/20180804-fig5f-data_AH2.pdf", sep=""), p3, width = 5, height = 6)


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
ggsave(paste("figure/20180804-fig5f-data_HN2.pdf", sep=""), p4, width = 5, height = 6)


