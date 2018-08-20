setwd("/mnt/bai/qinyuan/xiaoxuan/bac_New/")
print(paste("Your working directory is in",getwd()))
library(ggplot2)
library(dplyr)

load <- read.csv("table/20180802-fig5b.csv", header = T)

## AH-Bac: Dry & Wet, MH63 & WYJ, total microbiome load

t1 <- filter(load, Group == "AHBAC")

p1 <- ggplot(t1) +
  geom_boxplot(aes(x = Treatment, y = MicrobiomeLoad), width = 0.3, na.rm = TRUE) +
  geom_jitter(aes(x = Treatment, y = MicrobiomeLoad, color = Treatment), width = 0.08, na.rm = TRUE) +
  facet_wrap(~Genotype, ncol = 2) +
  labs(x="AHBAC-Group", y="Total microbiome load(relative to host plant )") +
  scale_colour_manual(values = c("goldenrod1", "deepskyblue")) +
  theme_bw()
p1
ggsave(paste("figure/20180802-totalload-ah-bac.pdf", sep=""), p1, width = 6, height = 5)


## HN-Bac: Dry & Wet, MH63 & WYJ, total microbiome load

t2 <- filter(load, Group == "HNBAC")

p2 <- ggplot(t2) +
  geom_boxplot(aes(x = Treatment, y = MicrobiomeLoad), width = 0.3, na.rm = TRUE) +
  geom_jitter(aes(x = Treatment, y = MicrobiomeLoad, color = Treatment), width = 0.08, na.rm = TRUE) +
  facet_wrap(~Genotype, ncol = 2) +
  labs(x="HNBAC-Group", y="Total microbiome load(relative to host plant )") +
  scale_colour_manual(values = c("goldenrod1", "deepskyblue")) +
  theme_bw()
p2
ggsave(paste("figure/20180802-totalload-hn-bac.pdf", sep=""), p2, width = 6, height = 5)


## AH-Fun: Dry & Wet, MH63 & WYJ, total microbiome load

t3 <- filter(load, Group == "AHFUN")

p3 <- ggplot(t3) +
  geom_boxplot(aes(x = Treatment, y = MicrobiomeLoad), width = 0.3, na.rm = TRUE) +
  geom_jitter(aes(x = Treatment, y = MicrobiomeLoad, color = Treatment), width = 0.08, na.rm = TRUE) +
  facet_wrap(~Genotype, ncol = 2) +
  labs(x="AHFUN-Group", y="Total microbiome load(relative to host plant )") +
  scale_colour_manual(values = c("goldenrod1", "deepskyblue")) +
  theme_bw()
p3
ggsave(paste("figure/20180802-totalload-ah-fun.pdf", sep=""), p3, width = 6, height = 5)


## HN-Bac: Dry & Wet, MH63 & WYJ, total microbiome load

t4 <- filter(load, Group == "HNFUN")

p4<- ggplot(t4) +
  geom_boxplot(aes(x = Treatment, y = MicrobiomeLoad), width = 0.3, na.rm = TRUE) +
  geom_jitter(aes(x = Treatment, y = MicrobiomeLoad, color = Treatment), width = 0.08, na.rm = TRUE) +
  facet_wrap(~Genotype, ncol = 2) +
  labs(x="HNFUN-Group", y="Total microbiome load(relative to host plant )") +
  scale_colour_manual(values = c("goldenrod1", "deepskyblue")) +
  theme_bw()
p4
ggsave(paste("figure/20180802-totalload-hn-fun.pdf", sep=""), p4, width = 6, height = 5)







