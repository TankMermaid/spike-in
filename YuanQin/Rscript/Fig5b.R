setwd("/mnt/bai/qinyuan/xiaoxuan/QA/")
print(paste("Your working directory is in",getwd()))
library(ggplot2)
library(dplyr)

load <- read.csv("data/Fig5b.csv", header = T)

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
ggsave(paste("result/Fig5b-AH-BAC.pdf", sep=""), p1, width = 6, height = 5)



## HN-Bac: Dry & Wet, MH63 & WYJ, total microbiome load
## Fig 5b

t2 <- filter(load, Group == "HNBAC")

p2 <- ggplot(t2) +
  geom_boxplot(aes(x = Treatment, y = MicrobiomeLoad), width = 0.3, na.rm = TRUE) +
  geom_jitter(aes(x = Treatment, y = MicrobiomeLoad, color = Treatment), width = 0.08, na.rm = TRUE) +
  facet_wrap(~Genotype, ncol = 2) +
  labs(x="HNBAC-Group", y="Total microbiome load(relative to host plant )") +
  scale_colour_manual(values = c("goldenrod1", "deepskyblue")) +
  theme_bw()
p2

ggsave(paste("result/Fig5b.pdf", sep=""), p2, width = 6, height = 5)


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
ggsave(paste("result/Fig5b-AH-Fun.pdf", sep=""), p3, width = 6, height = 5)


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
ggsave(paste("result/Fig5b-HN-Fun.pdf", sep=""), p4, width = 6, height = 5)





