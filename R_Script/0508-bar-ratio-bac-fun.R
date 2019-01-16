setwd("/mnt/bai/qinyuan/xiaoxuan/bac_New/")
print(paste("Your working directory is in",getwd()))
library(ggplot2)


## nature sample bacteria ratio aa & ra

bactable <- read.csv("table/0508ratio/20180508-ratio-bac-bar.csv", header = T)

bacplot <- ggplot(bactable) +
  geom_boxplot(aes(x = Approach, y = Ratio, fill = Approach), width = 0.3, na.rm = TRUE) + 
  geom_jitter(aes(x = Approach, y = Ratio, colour = Approach), width = 0.1, na.rm = TRUE) +
  scale_fill_hue(c = 30) +
  facet_wrap(~Group, ncol = 3) +
  scale_y_continuous(limits = c(1, 1.8)) + theme_bw() +
  labs(x="Approach", y="Ratio") 
bacplot


## nature sample fungi ratio aa & ra

funtable <- read.csv("table/0508ratio/20180508-ratio-fun-bar.csv", header = T)

funplot <- ggplot(funtable) +
  geom_boxplot(aes(x = Approach, y = Ratio, fill = Approach), width = 0.3, na.rm = TRUE) + 
  geom_jitter(aes(x = Approach, y = Ratio, colour = Approach), width = 0.1, na.rm = TRUE) +
  scale_fill_hue(c = 30) +
  facet_wrap(~Group, ncol = 3) +
  scale_y_continuous(limits = c(0, 22)) +
  theme_bw() +
  labs(x="Approach", y="Ratio") 
funplot
