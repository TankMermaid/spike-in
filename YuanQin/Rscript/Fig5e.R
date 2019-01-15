setwd("/mnt/bai/qinyuan/xiaoxuan/QA/")
print(paste("Your working directory is in",getwd()))

library(ggplot2)
library(dplyr)

data_ra <- read.csv("data/Fig5e-otu-bar-ra.csv", header = T)
data_aa <- read.csv("data/Fig5e-otu-bar-aa.csv", header = T)

## OTU11-RA
data1 <- filter(data_ra, OTUS == "OTU_11")

p1 <- ggplot(data1) +
  geom_boxplot(aes(x = Treatment, y = RA), width = 0.3, na.rm = TRUE) +
  geom_jitter(aes(x = Treatment, y = RA, color = Treatment), width = 0.08, na.rm = TRUE) +
  labs(x="Group_OTU_11", y="Relative genus abundance (RA) (‰)") +
  scale_colour_manual(values = c("goldenrod1", "deepskyblue")) +
  theme_bw()
p1
ggsave(paste("result/Fig5e-otus-11-ra.pdf", sep=""), p1, width = 4, height = 5)


## OTU16-RA
data2 <- filter(data_ra, OTUS == "OTU_16")

p2 <- ggplot(data2) +
  geom_boxplot(aes(x = Treatment, y = RA), width = 0.3, na.rm = TRUE) +
  geom_jitter(aes(x = Treatment, y = RA, color = Treatment), width = 0.08, na.rm = TRUE) +
  labs(x="Group_OTU_16", y="Relative genus abundance (RA) (‰)") +
  scale_colour_manual(values = c("goldenrod1", "deepskyblue")) +
  theme_bw()
p2
ggsave(paste("result/Fig5e-otus-16-ra.pdf", sep=""), p2, width = 4, height = 5)


## OTU13-RA
data3 <- filter(data_ra, OTUS == "OTU_13")

p3 <- ggplot(data3) +
  geom_boxplot(aes(x = Treatment, y = RA), width = 0.3, na.rm = TRUE) +
  geom_jitter(aes(x = Treatment, y = RA, color = Treatment), width = 0.08, na.rm = TRUE) +
  labs(x="Group_OTU_13", y="Relative genus abundance (RA) (‰)") +
  scale_colour_manual(values = c("goldenrod1", "deepskyblue")) +
  theme_bw()
p3
ggsave(paste("result/Fig5e-otus-13-ra.pdf", sep=""), p3, width = 4, height = 5)

## OTU11-AA

data4 <- filter(data_aa, OTUS == "OTU_11")

p4 <- ggplot(data4) +
  geom_boxplot(aes(x = Treatment, y = AA), width = 0.3, na.rm = TRUE) +
  geom_jitter(aes(x = Treatment, y = AA, color = Treatment), width = 0.08, na.rm = TRUE) +
  labs(x="Group_OTU_11", y="Quantitative genus abundance (QA)") +
  scale_colour_manual(values = c("goldenrod1", "deepskyblue")) +
  theme_bw()
p4
ggsave(paste("result/Fig5e-otus-11-aa.pdf", sep=""), p4, width = 4, height = 5)


## OTU16-AA

data5 <- filter(data_aa, OTUS == "OTU_16")

p5 <- ggplot(data5) +
  geom_boxplot(aes(x = Treatment, y = AA), width = 0.3, na.rm = TRUE) +
  geom_jitter(aes(x = Treatment, y = AA, color = Treatment), width = 0.08, na.rm = TRUE) +
  labs(x="Group_OTU_16", y="Quantitative genus abundance (QA)") +
  scale_colour_manual(values = c("goldenrod1", "deepskyblue")) +
  theme_bw()
p5
ggsave(paste("result/Fig5e-otus-16-aa.pdf", sep=""), p5, width = 4, height = 5)


## OTU13-AA

data6 <- filter(data_aa, OTUS == "OTU_13")

p6 <- ggplot(data5) +
  geom_boxplot(aes(x = Treatment, y = AA), width = 0.3, na.rm = TRUE) +
  geom_jitter(aes(x = Treatment, y = AA, color = Treatment), width = 0.08, na.rm = TRUE) +
  labs(x="Group_OTU_13", y="Quantitative genus abundance (QA)") +
  scale_colour_manual(values = c("goldenrod1", "deepskyblue")) +
  theme_bw()
p6
ggsave(paste("result/Fig5e-otus-13-aa.pdf", sep=""), p6, width = 4, height = 5)


## OTU12-RA
data7 <- filter(data_ra, OTUS == "OTU_12")

p7 <- ggplot(data7) +
  geom_boxplot(aes(x = Treatment, y = RA), width = 0.3, na.rm = TRUE) +
  geom_jitter(aes(x = Treatment, y = RA, color = Treatment), width = 0.08, na.rm = TRUE) +
  labs(x="Group_OTU_12", y="Relative genus abundance (RA) (‰)") +
  scale_colour_manual(values = c("goldenrod1", "deepskyblue")) +
  theme_bw()
p7
ggsave(paste("result/Fig5e-otus-12-ra.pdf", sep=""), p7, width = 4, height = 5)


## OTU20-RA
data8 <- filter(data_ra, OTUS == "OTU_20")

p8 <- ggplot(data8) +
  geom_boxplot(aes(x = Treatment, y = RA), width = 0.3, na.rm = TRUE) +
  geom_jitter(aes(x = Treatment, y = RA, color = Treatment), width = 0.08, na.rm = TRUE) +
  labs(x="Group_OTU_20", y="Relative genus abundance (RA) (‰)") +
  scale_colour_manual(values = c("goldenrod1", "deepskyblue")) +
  theme_bw()
p8
ggsave(paste("result/Fig5e-otus-20-ra.pdf", sep=""), p8, width = 4, height = 5)


## OTU12-AA

data9 <- filter(data_aa, OTUS == "OTU_12")

p9 <- ggplot(data9) +
  geom_boxplot(aes(x = Treatment, y = AA), width = 0.3, na.rm = TRUE) +
  geom_jitter(aes(x = Treatment, y = AA, color = Treatment), width = 0.08, na.rm = TRUE) +
  labs(x="Group_OTU_12", y="Quantitative genus abundance (QA)") +
  scale_colour_manual(values = c("goldenrod1", "deepskyblue")) +
  theme_bw()
p9
ggsave(paste("result/Fig5e-otus-12-aa.pdf", sep=""), p9, width = 4, height = 5)


## OTU20-AA

data10 <- filter(data_aa, OTUS == "OTU_20")

p10 <- ggplot(data10) +
  geom_boxplot(aes(x = Treatment, y = AA), width = 0.3, na.rm = TRUE) +
  geom_jitter(aes(x = Treatment, y = AA, color = Treatment), width = 0.08, na.rm = TRUE) +
  labs(x="Group_OTU_20", y="Quantitative genus abundance (QA)") +
  scale_colour_manual(values = c("goldenrod1", "deepskyblue")) +
  theme_bw()
p10
ggsave(paste("result/Fig5e-otus-20-aa.pdf", sep=""), p10, width = 4, height = 5)
