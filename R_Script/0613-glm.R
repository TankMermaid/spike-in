setwd("/mnt/bai/qinyuan/xiaoxuan/bac_New/")
print(paste("Your working directory is in",getwd()))
library(dplyr)
library(ggplot2)
library(broom)

spike <- read.csv("table/20180614-glm.csv", header = T)

# log10

mod1 <- glm(counts ~ log10, family = poisson(link = "log"), data = spike)

p1 <- ggplot(spike, aes(log10, counts, color = Dilution)) +
  geom_point(size = 2.5) +
  theme_bw() +
  scale_y_continuous(limits = c(0, 12000)) +
  geom_smooth(method = "glm", size = 0.8, colour = "grey10")+
  labs(x="spike-in amount (log10 copies/reaction)", y="spike-in read counts") 

p1


cor.test(spike$log10,spike$counts,method="pearson")

lm.fit1 <- lm(spike$counts~ spike$log10)
coef(lm.fit1)

ggsave(paste("figure/20180614-glm-log10.pdf", sep=""), p1, width = 7, height = 6)


# log2

mod2 <- glm(counts ~ log2, family = poisson(link = "log"), data = spike)

p2 <- ggplot(spike, aes(log2, counts, color = Dilution)) +
  geom_point(size = 2.5) +
  theme_bw() +
  scale_y_continuous(limits = c(0, 12000)) +
  geom_smooth(method = "glm", size = 0.8, colour = "grey10")+
  labs(x="spike-in amount (log2 copies/reaction)", y="spike-in read counts") 

p2



cor.test(spike$log2,spike$counts,method="pearson")

lm.fit1 <- lm(spike$counts~ spike$log2)
coef(lm.fit1)

s2 <- summary(mod2)


ggsave(paste("figure/20180614-glm-log2.pdf", sep=""), p2, width = 7, height = 6)
