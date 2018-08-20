setwd("/mnt/bai/qinyuan/xiaoxuan/QA/")
print(paste("Your working directory is in",getwd()))

library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyr)
library(dplyr)

fun_ra <- read.csv("data/Fig3-bar-ra.csv", header = T)
table_fun_ra_r <- gather(fun_ra,  key = sample , value = count, `sample01`:`sample05`)
table_fun_ra <- table_fun_ra_r
table_fun_ra$fungi <- factor(table_fun_ra_r$fungi)
table_fun_ra$group <- factor(table_fun_ra_r$group)
E0520_fun_ra <- filter(table_fun_ra, spike == "E0520")
group1_E0520_fun_ra <- filter(E0520_fun_ra , group == "1")
group2_E0520_fun_ra <- filter(E0520_fun_ra , group == "2")
group3_E0520_fun_ra <- filter(E0520_fun_ra , group == "3")

fun_aa <- read.csv("data/Fig3-bar-aa.csv", header = T)
table_fun_aa_r <- gather(fun_aa,  key = sample , value = count, `sample01`:`sample05`)
table_fun_aa <- table_fun_aa_r
table_fun_aa$fungi <- factor(table_fun_aa_r$fungi)
table_fun_aa$group <- factor(table_fun_aa_r$group)
E0520_fun_aa <- filter(table_fun_aa, spike == "E0520")
group1_E0520_fun_aa <- filter(E0520_fun_aa, group == "1")
group2_E0520_fun_aa <- filter(E0520_fun_aa, group == "2")
group3_E0520_fun_aa <- filter(E0520_fun_aa, group == "3")

Fig3c_ra_g1g2 <- rbind(group1_E0520_fun_ra, group2_E0520_fun_ra)
Fig3d_aa_g1g2 <- rbind(group1_E0520_fun_aa, group2_E0520_fun_aa)
Fig3e_ra_g1g3 <- rbind(group1_E0520_fun_ra, group3_E0520_fun_ra)
Fig3f_aa_g1g3 <- rbind(group1_E0520_fun_aa, group3_E0520_fun_aa)



Fig3c <- ggplot(Fig3c_ra_g1g2) +
  geom_boxplot(aes(x = fungi, y = count, fill = group), width = 0.3, na.rm = TRUE) +
  scale_fill_manual(values=c("#9acd32","#ffa07a"))+
  geom_line(aes(x = fungi, y = mean, group = group, color = group), na.rm = TRUE) +
  labs(x="Fungi", y="Relative abundance (%)")  +
  scale_x_discrete(limits=c("Basi-AF78","Asco-AF1","Asco-AF105"))+
  scale_color_manual(values=c("#9acd32","#ffa07a"))+
  theme_classic() +
  theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust = 0))
Fig3c

ggsave(paste("result/Fig3c-ra.pdf", sep=""), Fig3c, width = 3, height = 4)


Fig3d <- ggplot(Fig3d_aa_g1g2) +
  geom_boxplot(aes(x = fungi, y = count, fill = group), width = 0.3, na.rm = TRUE) +
  scale_fill_manual(values=c("#9acd32","#ffa07a"))+
  geom_line(aes(x = fungi, y = mean, group = group, color = group), na.rm = TRUE) +
  labs(x="Fungi", y="Quantitative abundance (relative to plant)")  +
  scale_x_discrete(limits=c("Basi-AF78","Asco-AF1","Asco-AF105"))+
  scale_color_manual(values=c("#9acd32","#ffa07a"))+
  theme_classic() +
  theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust = 0))
Fig3d

ggsave(paste("result/Fig3d-aa.pdf", sep=""), Fig3d, width = 3, height = 4)

Fig3e <- ggplot(Fig3e_ra_g1g3) +
  geom_boxplot(aes(x = fungi, y = count, fill = group), width = 0.3, na.rm = TRUE) +
  scale_fill_manual(values=c("#9acd32","#87cefa"))+
  geom_line(aes(x = fungi, y = mean, group = group, color = group), na.rm = TRUE) +
  labs(x="Fungi", y="Relative abundance (%)")  +
  scale_x_discrete(limits=c("Basi-AF78","Asco-AF1","Asco-AF105"))+
  scale_color_manual(values=c("#9acd32","#87cefa"))+
  theme_classic() +
  theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust = 0))
Fig3e

ggsave(paste("result/Fig3e-ra.pdf", sep=""), Fig3e, width = 3, height = 4)

Fig3f <- ggplot(Fig3f_aa_g1g3) +
  geom_boxplot(aes(x = fungi, y = count, fill = group), width = 0.3, na.rm = TRUE) +
  scale_fill_manual(values=c("#9acd32","#87cefa"))+
  geom_line(aes(x = fungi, y = mean, group = group, color = group), na.rm = TRUE) +
  labs(x="Fungi", y="Quantitative abundance (relative to plant)")  +
  scale_x_discrete(limits=c("Basi-AF78","Asco-AF1","Asco-AF105"))+
  scale_color_manual(values=c("#9acd32","#87cefa"))+
  theme_classic() +
  theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust = 0))
Fig3f

ggsave(paste("result/Fig3f-aa.pdf", sep=""), Fig3f, width = 3, height = 4)

