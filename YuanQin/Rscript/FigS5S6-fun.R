setwd("/mnt/bai/qinyuan/xiaoxuan/QA/")
print(paste("Your working directory is in",getwd()))
library(ggplot2)
library(tidyr)
library(dplyr)

ra<- read.csv("data/FigS5S6-ra-fun.csv", header = T)
aa<- read.csv("data/FigS5S6-aa-fun.csv", header = T)
group1 <- filter(ra, group == "1")
group2 <- filter(ra, group == "2")
group3 <- filter(ra, group == "3")

group1aa <- filter(aa, group == "1")
group2aa <- filter(aa, group == "2")
group3aa <- filter(aa, group == "3")

## fungi ra group1 & group2

data1 <- rbind(group1, group2)
colnames(data1)[1] = "fungi"
table1 <- gather(data1, key = sample , value = count, `sample01`:`sample05`)

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

ggsave(paste("result/FigS5-ra.pdf", sep=""), boxplot1, width = 9, height = 5)


## fungi ra group2-group3

data2 <- rbind(group2, group3)
colnames(data2)[1] = "fungi"
table2 <- gather(data2, key = sample , value = count, `sample01`:`sample05`)

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

ggsave(paste("result/FigS6-ra.pdf", sep=""), boxplot2, width = 9, height = 5)




## fungi aa group1-group2

data3 <- rbind(group1aa, group2aa)
colnames(data3)[1] = "fungi"
table3 <- gather(data3, key = sample , value = count, `sample01`:`sample05`)

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

ggsave(paste("result/FigS5-aa.pdf", sep=""), boxplot3, width = 9, height = 5)


## fungi spike aa group2-group3

data4 <- rbind(group3aa, group2aa)
colnames(data4)[1] = "fungi"
table4 <- gather(data4, key = sample , value = count, `sample01`:`sample05`)

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
ggsave(paste("result/FigS6-aa.pdf", sep=""), boxplot4, width = 9, height = 5)






