setwd("/mnt/bai/qinyuan/xiaoxuan/bac_New/")
print(paste("Your working directory is in",getwd()))
library(ggplot2)
library(tidyr)


## bacteria spike ra group2-group3

rawtable1 <- read.csv("table/20180504-spike4-ra-group2-3.csv", header = T)
colnames(rawtable1)[1] = "bacteria"
table1 <- gather(rawtable1, key = sample , value = count, `sample01`:`sample05`)

ratable1 <- table1
ratable1$bacteria <- factor(table1$bacteria)
ratable1$group <- factor(table1$group)

boxplot1 <- ggplot(ratable1) +
  geom_boxplot(aes(x = bacteria, y = count, fill = group), width = 0.3, na.rm = TRUE) +
  geom_line(aes(x = bacteria, y = mean, group = group, color = group), na.rm = TRUE) +
  facet_wrap(~spike, ncol = 4) +
  labs(x="Bacteria", y="Relative abundance (%)") +
  theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust = 0))
boxplot1

ggsave(paste("figure/20180504-spike4-ra-group2-3.pdf", sep=""), boxplot1, width = 9, height = 5)


## bacteria spike ra group1-group2

rawtable2 <- read.csv("table/20180504-bac-spike4-ra-group2-1.csv", header = T)
colnames(rawtable2)[1] = "bacteria"
table2 <- gather(rawtable2, key = sample , value = count, `sample01`:`sample05`)

ratable2 <- table2
ratable2$bacteria <- factor(table2$bacteria)
ratable2$group <- factor(table2$group)

boxplot2 <- ggplot(ratable2) +
  geom_boxplot(aes(x = bacteria, y = count, fill = group), width = 0.3, na.rm = TRUE) +
  geom_line(aes(x = bacteria, y = mean, group = group, color = group), na.rm = TRUE) +
  facet_wrap(~spike, ncol = 4) +
  labs(x="Bacteria", y="Relative abundance (%)") +
  theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust = 0))
boxplot2

ggsave(paste("figure/20180504-spike4-ra-group2-1.pdf", sep=""), boxplot2, width = 9, height = 5)


## bacteria spike aa group2-group3

rawtable3 <- read.csv("table/20180504-spike3-aa-group2-3.csv", header = T)
colnames(rawtable3)[1] = "bacteria"
table3 <- gather(rawtable3, key = sample , value = count, `sample01`:`sample05`)

aatable1 <- table3
aatable1$bacteria <- factor(table3$bacteria)
aatable1$group <- factor(table3$group)

boxplot3 <- ggplot(aatable1) +
  geom_boxplot(aes(x = bacteria, y = count, fill = group), width = 0.3, na.rm = TRUE) +
  geom_line(aes(x = bacteria, y = mean, group = group, color = group), na.rm = TRUE) +
  facet_wrap(~spike, ncol = 4) +
  labs(x="Bacteria", y="Absolute abundance") +
  theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust = 0))
boxplot3

ggsave(paste("figure/20180504-spike3-aa-group2-3.pdf", sep=""), boxplot1, width = 7, height = 5)


## bacteria spike aa group2-group1

rawtable4 <- read.csv("table/20180504-bac-spike3-aa-group2-1.csv", header = T)
colnames(rawtable4)[1] = "bacteria"
table4 <- gather(rawtable4, key = sample , value = count, `sample01`:`sample05`)

aatable2 <- table4
aatable2$bacteria <- factor(table4$bacteria)
aatable2$group <- factor(table4$group)

boxplot4 <- ggplot(aatable2) +
  geom_boxplot(aes(x = bacteria, y = count, fill = group), width = 0.3, na.rm = TRUE) +
  geom_line(aes(x = bacteria, y = mean, group = group, color = group), na.rm = TRUE) +
  facet_wrap(~spike, ncol = 4) +
  labs(x="Bacteria", y="Absolute abundance") +
  theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust = 0))
boxplot4

ggsave(paste("figure/20180504-spike3-aa-group2-1.pdf", sep=""), boxplot4, width = 7, height = 5)

