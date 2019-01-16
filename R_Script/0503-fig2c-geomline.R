setwd("/mnt/bai/qinyuan/xiaoxuan/bac_New/")
print(paste("Your working directory is in",getwd()))
library(ggplot2)
library(tidyr)

rawtable2 <- read.csv("table/fig2c/20180503-e0520-ra-group1-3.csv", header = T)
colnames(rawtable2)[1] = "bacteria"
rawtable4 <- read.csv("table/fig2c/20180503-e0520-aa-group1-3.csv", header = T)
colnames(rawtable4)[1] = "bacteria"

ra2 <- gather(rawtable2, key = sample , value = count, `sample01`:`sample05`)
aa2 <- gather(rawtable4, key = sample , value = count, `sample01`:`sample05`)

ra2data <- ra2
ra2data$bacteria <- factor(ra2data$bacteria)
ra2data$group <- factor(ra2data$group)
aa2data <- aa2
aa2data$bacteria <- factor(aa2data$bacteria)
aa2data$group <- factor(aa2data$group)
plot2 <- ggplot(ra2data) +
  geom_point(aes(x = bacteria, y = mean, color = group), na.rm = TRUE)+
  geom_line(aes(x = bacteria, y = mean, group = group, color = group), na.rm = TRUE) +
  labs(x="Bacteria", y="Relative abundance (%)") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust = 0))
plot2
plot4 <- ggplot(aa2data) +
  geom_point(aes(x = bacteria, y = mean, color = group), na.rm = TRUE)+
  geom_line(aes(x = bacteria, y = mean, group = group, color = group), na.rm = TRUE) +
  labs(x="Bacteria", y="Absolute abundance (%)") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust = 0))
plot4
ggsave(paste("figure/20180503-line-e0520-ra2.pdf", sep=""), plot2, width = 5, height = 5)
ggsave(paste("figure/20180503-line-e0520-aa2.pdf", sep=""), plot4, width = 5, height = 5)





