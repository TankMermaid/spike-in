setwd("/mnt/bai/qinyuan/xiaoxuan/QA/")
print(paste("Your working directory is in",getwd()))

library(ggplot2)
library(dplyr)

data_ra <- read.csv("data/Fig6d-ra.csv", header = T)
data_aa <- read.csv("data/Fig6d-aa.csv", header = T)


## RA

p1 <- ggplot(data_ra) +
  geom_boxplot(aes(x = Treatment, y = RA), width = 0.3, na.rm = TRUE) +
  geom_jitter(aes(x = Treatment, y = RA, color = Treatment), size = 2, width = 0.08, na.rm = TRUE) +
  labs(x="Group_OTU_11", y="Relative genus abundance (RA) (‰)") +
  scale_colour_manual(values = c("gray48", "olivedrab4"))+
  theme_classic()
p1
ggsave(paste("result/Fig6d-ra.pdf", sep=""), p1, width = 4, height = 5)


## AA

p6 <- ggplot(data_aa) +
  geom_boxplot(aes(x = Treatment, y = AA), width = 0.3, na.rm = TRUE) +
  geom_jitter(aes(x = Treatment, y = AA, color = Treatment), size = 2, width = 0.08, na.rm = TRUE) +
  labs(x="Group_OTU_13", y="Quantitative genus abundance (QA)") +
  scale_colour_manual(values = c("gray48", "olivedrab4")) +
  theme_classic()
p6
ggsave(paste("result/Fig6d-aa.pdf", sep=""), p6, width = 4, height = 5)








