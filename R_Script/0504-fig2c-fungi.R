setwd("/mnt/bai/qinyuan/xiaoxuan/bac_New/")
print(paste("Your working directory is in",getwd()))
library(ggplot2)
library(tidyr)

# fungi ra group 2-3 (1:1:1 & 2:2:1)

rawtable1<- read.csv("table/20180504-fun-ra-group2-3.csv", header = T)
colnames(rawtable1)[1] = "fungi"
table1 <- gather(rawtable1, key = sample , value = count, `sample01`:`sample05`)

ratable1 <- table1
ratable1$fungi <- factor(table1$fungi)
ratable1$group <- factor(table1$group)

raboxplot1 <- ggplot(ratable1) +
  geom_boxplot(aes(x = fungi, y = count, fill = group), width = 0.3, na.rm = TRUE) +
  geom_line(aes(x = fungi, y = mean, group = group, color = group), na.rm = TRUE) +
  facet_wrap(~spike, ncol = 5) +
  labs(x="Fungi", y="Relative abundance (%)") +
  theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust = 0))
raboxplot1

ggsave(paste("figure/20180504-fun-ra-group2-3.pdf", sep=""), raboxplot1, width = 9, height = 5)


# fungi ra group 1-2 (2:2:2 & 1:1:1)

rawtable2<- read.csv("table/20180504-fun-ra-group2-1.csv", header = T)
colnames(rawtable2)[1] = "fungi"
table2 <- gather(rawtable2, key = sample , value = count, `sample01`:`sample05`)

ratable2 <- table2
ratable2$fungi <- factor(table2$fungi)
ratable2$group <- factor(table2$group)

raboxplot2 <- ggplot(ratable2) +
  geom_boxplot(aes(x = fungi, y = count, fill = group), width = 0.3, na.rm = TRUE) +
  geom_line(aes(x = fungi, y = mean, group = group, color = group), na.rm = TRUE) +
  facet_wrap(~spike, ncol = 5) +
  labs(x="Fungi", y="Relative abundance (%)") +
  theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust = 0))
raboxplot2

ggsave(paste("figure/20180504-fun-ra-group2-1.pdf", sep=""), raboxplot1, width = 9, height = 5)

# fungi aa group 2-3 (1:1:1 & 2:2:1)

rawtable3<- read.csv("table/20180504-fun-aa-group2-3.csv", header = T)
colnames(rawtable3)[1] = "fungi"
table3 <- gather(rawtable3, key = sample , value = count, `sample01`:`sample05`)

aatable1 <- table3
aatable1$fungi <- factor(table3$fungi)
aatable1$group <- factor(table3$group)

aaboxplot1 <- ggplot(aatable1) +
  geom_boxplot(aes(x = fungi, y = count, fill = group), width = 0.3, na.rm = TRUE) +
  geom_line(aes(x = fungi, y = mean, group = group, color = group), na.rm = TRUE) +
  facet_wrap(~spike, ncol = 4) +
  labs(x="Fungi", y="Absolute abundance") +
  theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust = 0))
aaboxplot1

ggsave(paste("figure/20180504-fun-aa-group2-3.pdf", sep=""), aaboxplot1, width = 9, height = 5)


# fungi aa group 1-2 (2:2:2 & 2:2:1)

rawtable4<- read.csv("table/20180504-fun-aa-group2-1.csv", header = T)
colnames(rawtable4)[1] = "fungi"
table4 <- gather(rawtable4, key = sample , value = count, `sample01`:`sample05`)

aatable2 <- table4
aatable2$fungi <- factor(table4$fungi)
aatable2$group <- factor(table4$group)

aaboxplot2 <- ggplot(aatable2) +
  geom_boxplot(aes(x = fungi, y = count, fill = group), width = 0.3, na.rm = TRUE) +
  geom_line(aes(x = fungi, y = mean, group = group, color = group), na.rm = TRUE) +
  facet_wrap(~spike, ncol = 4) +
  labs(x="Fungi", y="Absolute abundance") +
  theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust = 0))
aaboxplot2

ggsave(paste("figure/20180504-fun-aa-group2-1.pdf", sep=""), aaboxplot2, width = 9, height = 5)





