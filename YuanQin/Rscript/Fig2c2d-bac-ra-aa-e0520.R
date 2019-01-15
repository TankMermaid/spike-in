setwd("/mnt/bai/qinyuan/xiaoxuan/QA/")
print(paste("Your working directory is in",getwd()))
library(ggplot2)
library(tidyr)
library(dplyr)

ra<- read.csv("data/Fig2c2d-ra-bac-e0520.csv", header = T)
aa<- read.csv("data/Fig2c2d-aa-bac-e0520.csv", header = T)
group1 <- filter(ra, group == "1")
group2 <- filter(ra, group == "2")
group3 <- filter(ra, group == "3")

group1aa <- filter(aa, group == "1")
group2aa <- filter(aa, group == "2")
group3aa <- filter(aa, group == "3")



## bacteria ra e05-20 group1 & group2

data1 <- rbind(group1, group2)
colnames(data1)[1] = "bacteria"
table1 <- gather(data1, key = sample , value = count, `sample01`:`sample05`)

ratable1 <- table1
ratable1$bacteria <- factor(table1$bacteria)
ratable1$group <- factor(table1$group)

boxplot1 <- ggplot(ratable1) +
  geom_boxplot(aes(x = bacteria, y = count, fill = group), width = 0.3, na.rm = TRUE) +
  geom_line(aes(x = bacteria, y = mean, group = group, color = group), na.rm = TRUE) +
  labs(x="Bacteria", y="Relative abundance (%)") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust = 0))
boxplot1

ggsave(paste("result/Fig2c-ra.pdf", sep=""), boxplot1, width = 6, height = 5)

## bacteria ra e05-20 group2 & group3

data2 <- rbind(group2, group3)
colnames(data2)[1] = "bacteria"
table2 <- gather(data2, key = sample , value = count, `sample01`:`sample05`)

ratable2 <- table2
ratable2$bacteria <- factor(table2$bacteria)
ratable2$group <- factor(table2$group)

boxplot2 <- ggplot(ratable2) +
  geom_boxplot(aes(x = bacteria, y = count, fill = group), width = 0.3, na.rm = TRUE) +
  geom_line(aes(x = bacteria, y = mean, group = group, color = group), na.rm = TRUE) +
  labs(x="Bacteria", y="Relative abundance (%)") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust = 0))
boxplot2

ggsave(paste("result/Fig2d-ra.pdf", sep=""), boxplot2, width = 6, height = 5)


## bacteria aa e05-20 group1 & group2

data3 <- rbind(group1aa, group2aa)
colnames(data3)[1] = "bacteria"
table3 <- gather(data3, key = sample , value = count, `sample01`:`sample05`)

aatable1 <- table3
aatable1$bacteria <- factor(table3$bacteria)
aatable1$group <- factor(table3$group)

boxplot3 <- ggplot(aatable1) +
  geom_boxplot(aes(x = bacteria, y = count, fill = group), width = 0.3, na.rm = TRUE) +
  geom_line(aes(x = bacteria, y = mean, group = group, color = group), na.rm = TRUE) +
  labs(x="Bacteria", y="Absolute abundance") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust = 0))
boxplot3

ggsave(paste("result/Fig2c-aa.pdf", sep=""), boxplot3, width = 6, height = 5)


## bacteria aa e05-20 group2 & group3

data4 <- rbind(group3aa, group2aa)
colnames(data4)[1] = "bacteria"
table4 <- gather(data4, key = sample , value = count, `sample01`:`sample05`)

aatable2 <- table4
aatable2$bacteria <- factor(table4$bacteria)
aatable2$group <- factor(table4$group)

boxplot4 <- ggplot(aatable2) +
  geom_boxplot(aes(x = bacteria, y = count, fill = group), width = 0.3, na.rm = TRUE) +
  geom_line(aes(x = bacteria, y = mean, group = group, color = group), na.rm = TRUE) +
  labs(x="Bacteria", y="Absolute abundance") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust = 0))
boxplot4
ggsave(paste("result/Fig2d-aa.pdf", sep=""), boxplot3, width = 6, height = 5)


