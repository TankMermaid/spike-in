setwd("/mnt/bai/qinyuan/xiaoxuan/bac_New/")
print(paste("Your working directory is in",getwd()))
library(ggplot2)
library(tidyr)

## fungi spike ra group1-group2

rawtable1 <- read.csv("table/20180619-fun-spike4-ra-group1-2.csv", header = T)
colnames(rawtable1)[1] = "fungi"
table1 <- gather(rawtable1, key = sample , value = count, `sample01`:`sample05`)

ratable1 <- table1
ratable1$fungi <- factor(table1$fungi)
ratable1$group <- factor(table1$group)

boxplot1 <- ggplot(ratable1) +
  geom_boxplot(aes(x = fungi, y = count, fill = group), width = 0.3, na.rm = TRUE) +
  geom_line(aes(x = fungi, y = mean, group = group, color = group), na.rm = TRUE) +
  facet_wrap(~spike, ncol = 3)  +
  scale_x_discrete(limits=c("Basi-AF78","Asco-AF1","Asco-AF105")) +
  labs(x="Fungi", y="Relative abundance (%)") + 
  theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust = 0))
boxplot1


ggsave(paste("figure/20180619-fun-spike4-ra-group1-2.pdf", sep=""), boxplot1, width = 9, height = 5)

## fungi spike ra group2-group3

rawtable2 <- read.csv("table/20180619-fun-spike4-ra-group2-3.csv", header = T)
colnames(rawtable2)[1] = "fungi"
table2 <- gather(rawtable2, key = sample , value = count, `sample01`:`sample05`)

ratable2 <- table2
ratable2$fungi <- factor(table2$fungi)
ratable2$group <- factor(table2$group)


boxplot2 <- ggplot(ratable2) +
  geom_boxplot(aes(x = fungi, y = count, fill = group), width = 0.3, na.rm = TRUE) +
  geom_line(aes(x = fungi, y = mean, group = group, color = group), na.rm = TRUE) +
  facet_wrap(~spike, ncol = 3)  +
  scale_x_discrete(limits=c("Basi-AF78","Asco-AF1","Asco-AF105")) +
  labs(x="Fungi", y="Relative abundance (%)") + 
  theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust = 0))
boxplot2


ggsave(paste("figure/20180619-fun-spike4-ra-group2-3.pdf", sep=""), boxplot2, width = 9, height = 5)


## fungi spike aa group1-group2

rawtable3 <- read.csv("table/20180619-fun-spike4-aa-group1-2.csv", header = T)
colnames(rawtable3)[1] = "fungi"
table3 <- gather(rawtable3, key = sample , value = count, `sample01`:`sample05`)

aatable1 <- table3
aatable1$fungi <- factor(table3$fungi)
aatable1$group <- factor(table3$group)

boxplot3 <- ggplot(aatable1) +
  geom_boxplot(aes(x = fungi, y = count, fill = group), width = 0.3, na.rm = TRUE) +
  geom_line(aes(x = fungi, y = mean, group = group, color = group), na.rm = TRUE) +
  facet_wrap(~spike, ncol = 3) +
  scale_x_discrete(limits=c("Basi-AF78","Asco-AF1","Asco-AF105"))+
  labs(x="Fungi", y="Quantitative abundance (relative to plant)") +
  theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust = 0))
boxplot3

ggsave(paste("figure/20180619-fun-spike4-aa-group1-2.pdf", sep=""), boxplot3, width = 9, height = 5)

## fungi spike aa group2-group3

rawtable4 <- read.csv("table/20180619-fun-spike4-aa-group2-3.csv", header = T)
colnames(rawtable4)[1] = "fungi"
table4 <- gather(rawtable4, key = sample , value = count, `sample01`:`sample05`)

aatable2 <- table4
aatable2$fungi <- factor(table4$fungi)
aatable2$group <- factor(table4$group)

boxplot4 <- ggplot(aatable2) +
  geom_boxplot(aes(x = fungi, y = count, fill = group), width = 0.3, na.rm = TRUE) +
  geom_line(aes(x = fungi, y = mean, group = group, color = group), na.rm = TRUE) +
  facet_wrap(~spike, ncol = 3) +
  scale_x_discrete(limits=c("Basi-AF78","Asco-AF1","Asco-AF105"))+
  labs(x="Fungi", y="Quantitative abundance (relative to plant)") +
  theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust = 0))
boxplot4
ggsave(paste("figure/20180619-fun-spike4-aa-group2-3.pdf", sep=""), boxplot4, width = 9, height = 5)
