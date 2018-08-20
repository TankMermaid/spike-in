setwd("/mnt/bai/qinyuan/xiaoxuan/QA/")
print(paste("Your working directory is in",getwd()))
library(dplyr)
library(ggplot2)
library(broom)

spike <- read.csv("data/Fig4-glm.csv", header = T)

# log2

mod2 <- glm(counts ~ log2, family = poisson(link = "log"), data = spike)

p <- ggplot(spike, aes(log2, counts, color = Dilution)) +
  geom_point(size = 2.5) +
  theme_classic()+
  geom_smooth(method = "glm", size = 0.8, colour = "grey10")+
  labs(x="spike-in amount (log2 copies/reaction)", y="spike-in read counts") 

p


cor.test(spike$log2,spike$counts,method="pearson")

lm.fit2 <- lm(spike$counts~ spike$log2)
coef(lm.fit2)

ggsave(paste("result/Fig4-glm-log2.pdf", sep=""), p, width = 7, height = 6)

