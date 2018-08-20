setwd("/mnt/bai/qinyuan/xiaoxuan/bac_New/")
print(paste("Your working directory is in",getwd()))
library(dplyr)
library(ggplot2)
library(broom)

spike <- read.csv("table/20180509-glm-ra.csv", header = T)

spike2 <- spike %>%
  mutate(
    lcopy = log10(spike_in_copynumber)
  )

plot <- ggplot(spike2, aes(lcopy, spike_in_RA)) +
  geom_point() +
  theme_bw() +
  geom_smooth(method = "glm", size = 1, colour = "lightcoral")

plot

mod <- glm(formula =spike_in_RA ~ lcopy, data = spike2, family = p())
coef(summary(mod))
