setwd("/mnt/bai/qinyuan/xiaoxuan/QA/")
print(paste("Your working directory is in",getwd()))
library(dplyr)
library(ggplot2)
library(tidyr)

data <- read.csv("./data/Fig6c-load-TAD-TAH.csv")

p = ggplot(data, aes(x=Group, y = Mean)) + 
  geom_bar(stat = "identity",fill = "white", colour = "black", width=0.7, na.rm = TRUE)+ 
  geom_errorbar(aes(ymin=Mean-Error, ymax=Mean+Error), width=.1) +
  xlab("Groups")+
  ylab("Quantitative abundance (relative to plant)")+ 
  theme_classic()
p
ggsave(paste("result/Fig6c-load-TAD-TAH.pdf", sep=""), p, width = 5, height =5)


