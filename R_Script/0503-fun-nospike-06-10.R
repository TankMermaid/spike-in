setwd("/mnt/bai/qinyuan/xiaoxuan/bac_New/")
print(paste("Your working directory is in",getwd()))
library(ggplot2)
library(tidyr)


data1 <- read.csv("table/20180503-fun-nospike-06-10.csv", header = T)
colnames(data1)[1] = "bacteria"

data2 <- gather(data1, key = replicate , value = count, `sample06`:`sample10`)

p <- ggplot(data2, aes(bacteria, count))+
  geom_boxplot(aes(fill = bacteria), na.rm = TRUE) +
  facet_wrap(~spike, ncol = 5) +
  theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust = 0)) +
  scale_fill_brewer(palette = "Set3") +
  labs(x="Bacteria", y="Count (%)")
p

ggsave(paste("figure/20180503-fun-nospike-06-10.pdf", sep=""), p, width = 5, height = 7)
