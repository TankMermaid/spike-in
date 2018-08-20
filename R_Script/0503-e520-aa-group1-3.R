setwd("/mnt/bai/qinyuan/xiaoxuan/bac_New/")
print(paste("Your working directory is in",getwd()))
library(ggplot2)
library(tidyr)


rawtable <- read.csv("table/fig2c/20180503-e0520-aa-group1-3.csv", header = T)
colnames(rawtable)[1] = "bacteria"
table <- gather(rawtable, key = sample , value = count, `sample01`:`sample05`)

table2 <- table
table2$bacteria <- factor(table2$bacteria)
table2$group <- factor(table2$group)

boxplot2 <- ggplot(table2) +
  geom_boxplot(aes(x = bacteria, y = count, fill = group), width = 0.3, na.rm = TRUE) +
  geom_line(aes(x = bacteria, y = mean, group = group, color = group), na.rm = TRUE) +
  labs(x="Bacteria", y="Relative abundance (%)") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust = 0))
boxplot2

ggsave(paste("figure/20180503-e0520-aa-group1-3.pdf", sep=""), boxplot2, width = 6, height = 5)
