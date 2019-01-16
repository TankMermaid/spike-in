setwd("/mnt/bai/qinyuan/xiaoxuan/bac_New/")
print(paste("Your working directory is in",getwd()))
library(dplyr)
library(ggplot2)
library(tidyr)


# bacteria aa E03

t1 <- read.csv("table/0606raaa/20180606-bacteria-aa-e03.csv", header = T)
bacaa03 <- gather(t1, key = variable , value = count, `one.fold`:`two.fold`)

bacaa03 = ggplot(bacaa03, aes(x=variable, y = count, fill = tax )) + 
  geom_bar(stat = "identity", width=0.7) + 
  scale_fill_manual(values = c("#323643", "#606470", "#93deff", "#bad7df","#ffe2e2", 
                               "#f6f6f6", "#99ddcc", "#ea4c88", "#eaafaf", "#a2738c", 
                               "#645c84", "#427996", "#fbe4c9", "#ff5d5d", "#952e4b", 
                               "#c06c84", "#f67280", "#f8b195","#ffe495","#ffc97b","#44558f"))+
  xlab("Groups")+
  ylab("Quantitative Abundance (relative to host plant)")+ 
  theme_classic()
bacaa03

ggsave(paste("figure/20180606-bacteria-aa-e03.pdf", sep=""), bacaa03, width = 6.5, height = 6)


# bacteria aa E04

t2 <- read.csv("table/0606raaa/20180606-bacteria-aa-e04.csv", header = T)
bacaa04 <- gather(t2, key = variable , value = count, `one.fold`:`two.fold`)

bacaa04 = ggplot(bacaa04, aes(x=variable, y = count, fill = tax )) + 
  geom_bar(stat = "identity", width=0.7) + 
  scale_fill_manual(values = c("#323643", "#606470", "#93deff", "#bad7df","#ffe2e2", 
                               "#f6f6f6", "#99ddcc", "#ea4c88", "#eaafaf", "#a2738c", 
                               "#645c84", "#427996", "#fbe4c9", "#ff5d5d", "#952e4b", 
                               "#c06c84", "#f67280", "#f8b195","#ffe495","#ffc97b","#44558f"))+
  xlab("Groups")+ylab("Quantitative Abundance (relative to host plant)")+ theme_classic()
bacaa04

ggsave(paste("figure/20180606-bacteria-aa-e04.pdf", sep=""), bacaa04, width = 6.5, height = 6)



# bacteria aa E55

t3 <- read.csv("table/0606raaa/20180606-bacteria-aa-e55.csv", header = T)
bacaa55 <- gather(t3, key = variable , value = count, `one.fold`:`two.fold`)

bacaa55 = ggplot(bacaa55, aes(x=variable, y = count, fill = tax )) + 
  geom_bar(stat = "identity", width=0.7) + 
  scale_fill_manual(values = c("#323643", "#606470", "#93deff", "#bad7df","#ffe2e2", 
                               "#f6f6f6", "#99ddcc", "#ea4c88", "#eaafaf", "#a2738c", 
                               "#645c84", "#427996", "#fbe4c9", "#ff5d5d", "#952e4b", 
                               "#c06c84", "#f67280", "#f8b195","#ffe495","#ffc97b","#44558f"))+
  xlab("Groups")+ylab("Quantitative Abundance (relative to host plant)")+ theme_classic()
bacaa55

ggsave(paste("figure/20180606-bacteria-aa-e55.pdf", sep=""), bacaa55, width = 6.5, height = 6)


# bacteria ra E03

t4 <- read.csv("table/0606raaa/20180606-bacteria-ra-e03.csv", header = T)
bacra03 <- gather(t4, key = variable , value = count, `one.fold`:`two.fold`)

bacra03 = ggplot(bacra03, aes(x=variable, y = count, fill = tax )) + 
  geom_bar(stat = "identity",position="fill", width=0.7)+ 
  scale_fill_manual(values = c("#323643", "#606470", "#93deff", "#bad7df","#ffe2e2", 
                               "#f6f6f6", "#99ddcc", "#ea4c88", "#eaafaf", "#a2738c", 
                               "#645c84", "#427996", "#fbe4c9", "#ff5d5d", "#952e4b", 
                               "#c06c84", "#f67280", "#f8b195","#ffe495","#ffc97b","#44558f"))+
  scale_y_continuous(labels = scales::percent) + 
  xlab("Groups")+ylab("Relative Abundance (%)")+ theme_classic()
bacra03

ggsave(paste("figure/20180606-bacteria-ra-e03.pdf", sep=""), bacra03, width =  6.5, height = 6)



# bacteria ra E04

t5 <- read.csv("table/0606raaa/20180606-bacteria-ra-e04.csv", header = T)
bacra04 <- gather(t5, key = variable , value = count, `one.fold`:`two.fold`)

bacra04 = ggplot(bacra04, aes(x=variable, y = count, fill = tax )) + 
  geom_bar(stat = "identity",position="fill", width=0.7)+ 
  scale_fill_manual(values = c("#323643", "#606470", "#93deff", "#bad7df","#ffe2e2", 
                               "#f6f6f6", "#99ddcc", "#ea4c88", "#eaafaf", "#a2738c", 
                               "#645c84", "#427996", "#fbe4c9", "#ff5d5d", "#952e4b", 
                               "#c06c84", "#f67280", "#f8b195","#ffe495","#ffc97b","#44558f"))+
  scale_y_continuous(labels = scales::percent) + 
  xlab("Groups")+ylab("Relative Abundance (%)")+ theme_classic()
bacra04

ggsave(paste("figure/20180606-bacteria-ra-e04.pdf", sep=""), bacra04, width = 6.5, height = 6)


# bacteria ra E55

t6 <- read.csv("table/0606raaa/20180606-bacteria-ra-e55.csv", header = T)
bacra55 <- gather(t6, key = variable , value = count, `one.fold`:`two.fold`)

bacra55 = ggplot(bacra55, aes(x=variable, y = count, fill = tax )) + 
  geom_bar(stat = "identity",position="fill", width=0.7)+ 
  scale_fill_manual(values = c("#323643", "#606470", "#93deff", "#bad7df","#ffe2e2", 
                               "#f6f6f6", "#99ddcc", "#ea4c88", "#eaafaf", "#a2738c", 
                               "#645c84", "#427996", "#fbe4c9", "#ff5d5d", "#952e4b", 
                               "#c06c84", "#f67280", "#f8b195","#ffe495","#ffc97b","#44558f"))+
  scale_y_continuous(labels = scales::percent) + 
  xlab("Groups")+ylab("Relative Abundance (%)")+ theme_classic()
bacra55

ggsave(paste("figure/20180606-bacteria-ra-e55.pdf", sep=""), bacra55, width = 6.5, height = 6)


# fungi aa E03

t7 <- read.csv("table/0606raaa/20180606-fungi-aa-e03.csv", header = T)
funaa03 <- gather(t7, key = variable , value = count, `one.fold`:`two.fold`)

funaa03 = ggplot(funaa03, aes(x=variable, y = count, fill = tax )) + 
  geom_bar(stat = "identity", width=0.7) + 
  scale_fill_manual(values = c("#323643", "#606470", "#93deff", "#bad7df","#ffe2e2", 
                               "#f6f6f6", "#99ddcc", "#ea4c88", "#eaafaf", "#a2738c", 
                               "#645c84"))+
  xlab("Groups")+ylab("Quantitative Abundance (relative to host plant)")+ 
  theme_classic()
funaa03

ggsave(paste("figure/20180606-fungi-aa-e03.pdf", sep=""), funaa03, width = 5, height = 6)


# fungi ra E03

t8 <- read.csv("table/0606raaa/20180606-fungi-ra-e03.csv", header = T)
funra03 <- gather(t8, key = variable , value = count, `one.fold`:`two.fold`)

funra03 = ggplot(funra03, aes(x=variable, y = count, fill = tax )) + 
  geom_bar(stat = "identity",position="fill", width=0.7)+ 
  scale_fill_manual(values = c("#323643", "#606470", "#93deff", "#bad7df","#ffe2e2", 
                               "#f6f6f6", "#99ddcc", "#ea4c88", "#eaafaf", "#a2738c", 
                               "#645c84"))+
  scale_y_continuous(labels = scales::percent) + 
  xlab("Groups")+ylab("Relative Abundance (%)")+ 
  theme_classic()
funra03

ggsave(paste("figure/20180606-fungi-ra-e03.pdf", sep=""), funra03, width = 5, height = 6)




