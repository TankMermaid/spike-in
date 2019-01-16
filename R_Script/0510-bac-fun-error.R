setwd("/mnt/bai/qinyuan/xiaoxuan/bac_New/")
print(paste("Your working directory is in",getwd()))
library(ggplot2)


## mock bacteria error

bac_error <- read.csv("table/20180510-bac-error.csv", header = T)

bacplot <- ggplot(bac_error) +
  geom_boxplot(aes(x = approach, y = error, fill = approach), width = 0.3, na.rm = TRUE) + 
  geom_jitter(aes(x = approach, y = error, colour = approach), width = 0.1, na.rm = TRUE) +
  scale_fill_hue(c = 30) +
  facet_wrap(~spike, ncol = 4) +
  theme_bw() +
  labs(x="Approach", y="Error") 
bacplot

ggsave(paste("figure/2018010-bac-error.pdf", sep=""), bacplot, width = 7, height = 5)





## mock fungi error

fun_error <- read.csv("table/20180510-fun-error.csv", header = T)

funplot <- ggplot(fun_error) +
  geom_boxplot(aes(x = approach, y = error, fill = approach), width = 0.3, na.rm = TRUE) + 
  geom_jitter(aes(x = approach, y = error, colour = approach), width = 0.1, na.rm = TRUE) +
  scale_fill_hue(c = 30) +
  facet_wrap(~spike, ncol = 4) +
  theme_bw() +
  labs(x="Approach", y="Error") 
funplot

ggsave(paste("figure/2018010-fun-error.pdf", sep=""), funplot, width = 7, height = 5)









