setwd("/mnt/bai/qinyuan/xiaoxuan/bac_New/")
print(paste("Your working directory is in",getwd()))
library(ggplot2)
library(tidyr)

spike <- read.csv("table/20180413-bac-spike.csv", header = T)
colnames(spike)[1] = "sample"
spikedata <- data.frame(spike$sample, spike$E05, spike$E05.5, spike$E05.10, spike$E05.20, spike$E00, spike$group)
colnames(spikedata)[2] = "E05"
colnames(spikedata)[3] = "E05-5"
colnames(spikedata)[4] = "E05-10"
colnames(spikedata)[5] = "E05-20"
colnames(spikedata)[6] = "E00"
colnames(spikedata)[7] = "Group"

data <- gather(spikedata, key = spikein , value = count, `E05`:`E00`)


p <- ggplot(data, aes(spikein, count)) + 
  geom_boxplot(aes(fill = spikein), na.rm = TRUE) +
  facet_wrap(~Group, ncol = 3) +
  theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust = 0)) +
  scale_fill_brewer(palette = "Set3") +
  labs(x="spikein", y="Count (%)") +
  scale_x_discrete(limits=c("E00","E05","E05-5","E05-10","E05-20"))
p
ggsave(paste("figure/20180413-bac-spike.pdf", sep=""), p, width = 8, height = 5)
