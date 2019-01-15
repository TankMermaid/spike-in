setwd("/mnt/bai/qinyuan/xiaoxuan/QA/")
print(paste("Your working directory is in",getwd()))
library(ggplot2)
library(tidyr)

data <- read.csv("data/FigS4-fun-nospike-group1-group3.csv", header = T)
colnames(data)[1] = "fungi"
group1 <- select(data, fungi, spike, sample06, sample09, sample10)
group3 <- select(data, fungi, spike, sample11, sample12, sample13, sample14, sample15)


# fungi 1:1:1 (group1)

data1 <- gather(group1, key = replicate , value = count, `sample06`:`sample10`)

p1 <- ggplot(data1, aes(fungi, count))+
  geom_boxplot(aes(fill = fungi), na.rm = TRUE) +
  facet_wrap(~spike, ncol = 5)  +
  labs(x="Dilution gradient spike-in", y="Relative abundance (%)")+
  scale_fill_manual(values = c("#8ea6b4","#e7eff3", "#f8d5f0")) +
  scale_x_discrete(limits=c("Basi-AF78","Asco-AF1","Asco-AF105"))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust = 0)) 
p1

ggsave(paste("result/FigS4-fun-group1.pdf", sep=""), p1, width = 8, height = 5)

# fungi 2:2:1 (group3)

data3 <- gather(group3, key = replicate , value = count, `sample11`:`sample15`)


p2<- ggplot(data3, aes(fungi, count))+
  geom_boxplot(aes(fill = fungi), na.rm = TRUE) +
  facet_wrap(~spike, ncol = 5)  +
  labs(x="Dilution gradient spike-in", y="Relative abundance (%)")+
  scale_fill_manual(values = c("#8ea6b4","#e7eff3", "#f8d5f0")) +
  scale_x_discrete(limits=c("Basi-AF78","Asco-AF1","Asco-AF105"))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust = 0)) 
p2

ggsave(paste("result/FigS4-fun-group3.pdf", sep=""), p2, width = 8, height = 5)