setwd("/mnt/bai/qinyuan/xiaoxuan/bac_New/")
print(paste("Your working directory is in",getwd()))

library(ggplot2)
library(dplyr)

data <- read.csv("table/20180731-tahd-load.csv", header = T)



p <- ggplot(data) +
  geom_boxplot(aes(x = micro, y = load), width = 0.3, na.rm = TRUE) +
  geom_jitter(aes(x = micro, y = load, color = micro), width = 0.08, na.rm = TRUE) +
    labs( y="Quantitative abundance (relative toplant)") +
  scale_colour_manual(values = c("gray48", "olivedrab4", "royalblue4","plum")) +
  theme_bw()
p

data2 <- read.csv("table/20180731-tahd-load2.csv", header = T)



p2 <- ggplot(data2) +
  geom_boxplot(aes(x = micro, y = load), width = 0.3, na.rm = TRUE) +
  geom_jitter(aes(x = micro, y = load, color = micro), width = 0.08, na.rm = TRUE) +
  labs( y="Quantitative abundance (relative toplant)") +
  scale_colour_manual(values = c("gray48", "olivedrab4", "royalblue4","plum")) +
  theme_bw()
p2