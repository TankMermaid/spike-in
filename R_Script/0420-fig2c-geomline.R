setwd("/mnt/bai/qinyuan/xiaoxuan/bac_New/")
print(paste("Your working directory is in",getwd()))
library(ggplot2)
library(tidyr)


rawtable1 <- read.csv("table/20180418-e0520-group1-2.csv", header = T)
colnames(rawtable1)[1] = "bacteria"
rawtable2 <- read.csv("table/20180418-e0520-group1-3.csv", header = T)
colnames(rawtable2)[1] = "bacteria"
rawtable3 <- read.csv("table/20180420-e0520-group1-2.csv", header = T)
colnames(rawtable3)[1] = "bacteria"
rawtable4 <- read.csv("table/20180420-e0520-group1-3.csv", header = T)
colnames(rawtable4)[1] = "bacteria"

ra1 <- gather(rawtable1, key = sample , value = count, `sample01`:`sample05`)
ra2 <- gather(rawtable2, key = sample , value = count, `sample01`:`sample05`)
aa1 <- gather(rawtable3, key = sample , value = count, `sample01`:`sample05`)
aa2 <- gather(rawtable4, key = sample , value = count, `sample01`:`sample05`)


ra1data <- ra1
ra1data$bacteria <- factor(ra1data$bacteria)
ra1data$group <- factor(ra1data$group)

ra2data <- ra2
ra2data$bacteria <- factor(ra2data$bacteria)
ra2data$group <- factor(ra2data$group)

aa1data <- aa1
aa1data$bacteria <- factor(aa1data$bacteria)
aa1data$group <- factor(aa1data$group)

aa2data <- aa2
aa2data$bacteria <- factor(aa2data$bacteria)
aa2data$group <- factor(aa2data$group)



plot1 <- ggplot(ra1data) +
  geom_point(aes(x = bacteria, y = mean, color = group), na.rm = TRUE)+
  geom_line(aes(x = bacteria, y = mean, group = group, color = group), na.rm = TRUE) +
  labs(x="Bacteria", y="Relative abundance (%)") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust = 0))
plot1

plot2 <- ggplot(ra2data) +
  geom_point(aes(x = bacteria, y = mean, color = group), na.rm = TRUE)+
  geom_line(aes(x = bacteria, y = mean, group = group, color = group), na.rm = TRUE) +
  labs(x="Bacteria", y="Relative abundance (%)") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust = 0))
plot2


plot3 <- ggplot(aa1data) +
  geom_point(aes(x = bacteria, y = mean, color = group), na.rm = TRUE)+
  geom_line(aes(x = bacteria, y = mean, group = group, color = group), na.rm = TRUE) +
  labs(x="Bacteria", y="Absolute abundance (%)") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust = 0))
plot3

plot4 <- ggplot(aa2data) +
  geom_point(aes(x = bacteria, y = mean, color = group), na.rm = TRUE)+
  geom_line(aes(x = bacteria, y = mean, group = group, color = group), na.rm = TRUE) +
  labs(x="Bacteria", y="Absolute abundance (%)") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust = 0))
plot4


ggsave(paste("figure/20180420-line-e0520-ra1.pdf", sep=""), plot1, width = 5, height = 5)
ggsave(paste("figure/20180420-line-e0520-ra2.pdf", sep=""), plot2, width = 5, height = 5)
ggsave(paste("figure/20180420-line-e0520-aa1.pdf", sep=""), plot3, width = 5, height = 5)
ggsave(paste("figure/20180420-line-e0520-aa2.pdf", sep=""), plot4, width = 5, height = 5)
