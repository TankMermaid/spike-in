setwd("/mnt/bai/qinyuan/xiaoxuan/bac_New/")
print(paste("Your working directory is in",getwd()))

library(ggplot2)
library(dplyr)

data_ra <- read.csv("table/20187023-fig5d-ra-g.csv", header = T)
data_aa <- read.csv("table/20180723-fig5d-aa-g.csv", header = T)

## Intrasporangium_ra

data1 <- filter(data_ra, Genus == "Intrasporangium")

p1 <- ggplot(data1) +
  geom_boxplot(aes(x = Treatment, y = RA), width = 0.3, na.rm = TRUE) +
  geom_jitter(aes(x = Treatment, y = RA, color = Treatment), width = 0.08, na.rm = TRUE) +
  labs(x="Group_Intrasporangium", y="Relative genus abundance (RMP) (‰)") +
  scale_colour_manual(values = c("goldenrod1", "deepskyblue")) +
  theme_bw()
p1
ggsave(paste("figure/20180723-genus-in-ra.pdf", sep=""), p1, width = 4, height = 5)

## Intrasporangium-aa

data2 <- filter(data_aa, Genus == "Intrasporangium")

p2 <- ggplot(data2) +
  geom_boxplot(aes(x = Treatment, y = AA), width = 0.3, na.rm = TRUE) +
  geom_jitter(aes(x = Treatment, y = AA, color = Treatment), width = 0.08, na.rm = TRUE) +
  labs(x="Group_Intrasporangium", y="Quantitative genus abundance (QMP)") +
  scale_colour_manual(values = c("goldenrod1", "deepskyblue")) +
  theme_bw()
p2
ggsave(paste("figure/20180723-genus-in-aa.pdf", sep=""), p2, width = 4, height = 5)

## Microbacterium-ra

data3 <- filter(data_ra, Genus == "Microbacterium")

p3 <- ggplot(data3) +
  geom_boxplot(aes(x = Treatment, y = RA), width = 0.3, na.rm = TRUE) +
  geom_jitter(aes(x = Treatment, y = RA, color = Treatment), width = 0.08, na.rm = TRUE) +
  labs(x="Group_Microbacterium", y="Relative genus abundance (RMP) (‰)") +
  scale_colour_manual(values = c("goldenrod1", "deepskyblue")) +
  theme_bw()
p3
ggsave(paste("figure/20180723-genus-mi-ra.pdf", sep=""), p3, width = 4, height = 5)

## Microbacterium-aa

data4 <- filter(data_aa, Genus == "Microbacterium")

p4 <- ggplot(data4) +
  geom_boxplot(aes(x = Treatment, y = AA), width = 0.3, na.rm = TRUE) +
  geom_jitter(aes(x = Treatment, y = AA, color = Treatment), width = 0.08, na.rm = TRUE) +
  labs(x="Group_Microbacterium", y="Quantitative genus abundance (QMP)") +
  scale_colour_manual(values = c("goldenrod1", "deepskyblue")) +
  theme_bw()
p4
ggsave(paste("figure/20180723-genus-mi-aa.pdf", sep=""), p4, width = 4, height = 5)

## Mesorhizobium-ra

data5 <- filter(data_ra, Genus == "Mesorhizobium")

p5 <- ggplot(data5) +
  geom_boxplot(aes(x = Treatment, y = RA), width = 0.3, na.rm = TRUE) +
  geom_jitter(aes(x = Treatment, y = RA, color = Treatment), width = 0.08, na.rm = TRUE) +
  labs(x="Group_Mesorhizobium", y="Relative genus abundance (RMP) (‰)") +
  scale_colour_manual(values = c("goldenrod1", "deepskyblue")) +
  theme_bw()
p5
ggsave(paste("figure/20180723genus-me-ra.pdf", sep=""), p5, width = 4, height = 5)

## Mesorhizobium-aa

data6 <- filter(data_aa, Genus == "Mesorhizobium")

p6 <- ggplot(data6) +
  geom_boxplot(aes(x = Treatment, y = AA), width = 0.3, na.rm = TRUE) +
  geom_jitter(aes(x = Treatment, y = AA, color = Treatment), width = 0.08, na.rm = TRUE) +
  labs(x="Group_Mesorhizobium", y="Quantitative genus abundance (QMP)") +
  scale_colour_manual(values = c("goldenrod1", "deepskyblue")) +
  theme_bw()
p6
ggsave(paste("figure/20180723-genus-me-aa.pdf", sep=""), p6, width = 4, height = 5)

## Pelobacter_ra

data7 <- filter(data_ra, Genus == "Pelobacter")

p7 <- ggplot(data7) +
  geom_boxplot(aes(x = Treatment, y = RA), width = 0.3, na.rm = TRUE) +
  geom_jitter(aes(x = Treatment, y = RA, color = Treatment), width = 0.08, na.rm = TRUE) +
  labs(x="Group_Pelobacter", y="Relative genus abundance (RMP) (‰)") +
  scale_colour_manual(values = c("goldenrod1", "deepskyblue")) +
  theme_bw()
p7
ggsave(paste("figure/20180723genus-p-ra.pdf", sep=""), p7, width = 4, height = 5)

## Pelobacter-aa

data8 <- filter(data_aa, Genus == "Pelobacter")

p8 <- ggplot(data8) +
  geom_boxplot(aes(x = Treatment, y = AA), width = 0.3, na.rm = TRUE) +
  geom_jitter(aes(x = Treatment, y = AA, color = Treatment), width = 0.08, na.rm = TRUE) +
  labs(x="Group_Pelobacter", y="Quantitative genus abundance (QMP)") +
  scale_colour_manual(values = c("goldenrod1", "deepskyblue")) +
  theme_bw()
p8
ggsave(paste("figure/20180723-genus-p-aa.pdf", sep=""), p8, width = 4, height = 5)
