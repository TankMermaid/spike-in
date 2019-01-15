setwd("/mnt/bai/qinyuan/xiaoxuan/QA/")
print(paste("Your working directory is in",getwd()))
library(ggplot2)
library(tidyr)

# fungi 2:2:2

data1 <- read.csv("data/Fig3b-fun-nospike-group2.csv", header = T)
colnames(data1)[1] = "fungi"

data2 <- gather(data1, key = replicate , value = count, `sample01`:`sample05`)


p <- ggplot(data2, aes(fungi, count))+
  geom_boxplot(aes(fill = fungi), na.rm = TRUE) +
  facet_wrap(~spike, ncol = 5) +
  labs(x="Dilution gradient spike-in", y="Relative abundance (%)")+
  scale_fill_manual(values = c("#8ea6b4","#e7eff3", "#f8d5f0")) +
  scale_x_discrete(limits=c("Basi-AF78","Asco-AF1","Asco-AF105"))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust = 0)) 

p

ggsave(paste("result/Fig3b-fun-nospike-group2.pdf", sep=""), p, width = 8, height = 5)

