setwd("/mnt/bai/qinyuan/xiaoxuan/bac_New/")
print(paste("Your working directory is in",getwd()))
library(ggplot2)


## nature sample bacteria ratio aa & ra

table1 <- read.csv("table/0508ratio/20180508-ratio-bac.csv", header = T)

bactable <- table1
bactable$OTU.ID <- factor(table1$OTU.ID)
bactable$Group <- factor(table1$Group)
bactable$Approach <- factor(table1$Approach)

bacplot <- ggplot(bactable) +
  geom_point(aes(x = OTU.ID, y = Ratio, shape = Approach, color = Group), na.rm = TRUE)+
  geom_line(aes(x = OTU.ID, y = Ratio, group = Group, color = Group), na.rm = TRUE) +
  labs(x="OTU.ID", y="Ratio") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust = 0)) + 
  scale_y_continuous(limits = c(0, 2))
bacplot




## nature sample fungi ratio aa & ra

table2 <- read.csv("table/0508ratio/20180508-ratio-fun.csv", header = T)

funtable <- table2
funtable$OTU.ID <- factor(table2$OTU.ID)
funtable$Group <- factor(table2$Group)
funtable$Approach <- factor(table2$Approach)

funplot <- ggplot(funtable) +
  geom_point(aes(x = OTU.ID, y = Ratio, shape = Approach, color = Group), na.rm = TRUE)+
  geom_line(aes(x = OTU.ID, y = Ratio, group = Group, color = Group), na.rm = TRUE) +
  labs(x="OTU.ID", y="Ratio") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust = 0)) + 
  scale_y_continuous(limits = c(0, 30))
funplot
