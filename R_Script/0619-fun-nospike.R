setwd("/mnt/bai/qinyuan/xiaoxuan/bac_New/")
print(paste("Your working directory is in",getwd()))
library(ggplot2)
library(tidyr)

# fungi 2:2:2

data1 <- read.csv("table/20180619-fun-nospike-01-05.csv", header = T)
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

ggsave(paste("figure/20180619-fun-nospike-01-05.pdf", sep=""), p, width = 8, height = 5)


# fungi 1:1:1 

data3 <- read.csv("table/20180619-fun-nospike-06-10.csv", header = T)
colnames(data3)[1] = "fungi"

data4 <- gather(data3, key = replicate , value = count, `sample01`:`sample04`)


p2 <- ggplot(data4, aes(fungi, count))+
  geom_boxplot(aes(fill = fungi), na.rm = TRUE) +
  facet_wrap(~spike, ncol = 5)  +
  labs(x="Dilution gradient spike-in", y="Relative abundance (%)")+
  scale_fill_manual(values = c("#8ea6b4","#e7eff3", "#f8d5f0")) +
  scale_x_discrete(limits=c("Basi-AF78","Asco-AF1","Asco-AF105"))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust = 0)) 
p2

ggsave(paste("figure/20180619-fun-nospike-06-10.pdf", sep=""), p2, width = 8, height = 5)

# fungi 2:2:1 

data5 <- read.csv("table/20180619-fun-nospike-11-15.csv", header = T)
colnames(data5)[1] = "fungi"

data6 <- gather(data5, key = replicate , value = count, `sample01`:`sample05`)


p3 <- ggplot(data6, aes(fungi, count))+
  geom_boxplot(aes(fill = fungi), na.rm = TRUE) +
  facet_wrap(~spike, ncol = 5)  +
  labs(x="Dilution gradient spike-in", y="Relative abundance (%)")+
  scale_fill_manual(values = c("#8ea6b4","#e7eff3", "#f8d5f0")) +
  scale_x_discrete(limits=c("Basi-AF78","Asco-AF1","Asco-AF105"))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust = 0)) 
p3

ggsave(paste("figure/20180619-fun-nospike-11-15.pdf", sep=""), p3, width = 8, height = 5)

