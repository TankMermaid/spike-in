setwd("/mnt/bai/qinyuan/xiaoxuan/bac_New/")
print(paste("Your working directory is in",getwd()))
library(dplyr)
library(ggplot2)

## bac ra_aa error

spike <- read.csv("table/20180516-bac-error.csv", header = T)
x <- c(0, 1, 1.3)
y <- c(0, 1, 1.3)
df <- data.frame(x=x, y=y)

ggplot(df, aes(x, y)) +
  geom_line()


bacplot <- ggplot(df, aes(x, y)) +
  geom_line() +
  geom_point(data = spike, aes(x = error_RA, y = error_AA, color = Group.2), na.rm = TRUE) + 
  theme_bw() +  
  scale_x_continuous(limits = c(0, 1.3)) + 
  scale_y_continuous(limits = c(0, 1.3)) + 
  labs(x="RA error", y="AA error") 
bacplot


## fun ra_aa error

spike2 <- read.csv("table/20180619-fun-error.csv", header = T)

funplot <- ggplot(df, aes(x, y)) +
  geom_line() +
  geom_point(data = spike2, aes(x = RAerror, y = AAerror, color = Group.2), na.rm = TRUE) + 
  theme_bw() +  
  scale_x_continuous(limits = c(0, 1.3)) + 
  scale_y_continuous(limits = c(0, 1.3)) + 
  labs(x="Fungi RA error", y="Fungi AA error") 
funplot









