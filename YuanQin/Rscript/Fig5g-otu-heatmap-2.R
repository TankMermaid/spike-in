setwd("/mnt/bai/qinyuan/xiaoxuan/QA/")
print(paste("Your working directory is in",getwd()))

library(ggplot2)
library(plyr)
library(scales)
library(reshape2)
library(pheatmap)

data_1 <- read.csv("data/Fig5g-otu-heatmap-2-hn-MH-RA.csv", header = T)
data_2 <- read.csv("data/Fig5g-otu-heatmap-2-hn-MH-AA.csv", header = T)
data_3 <- read.csv("data/Fig5g-otu-heatmap-2-hn-WYJ-RA.csv", header = T)
data_4 <- read.csv("data/Fig5g-otu-heatmap-2-hn-WYJ-AA.csv", header = T)
data_5 <- read.csv("data/Fig5g-otu-heatmap-2-ah-MH-RA.csv", header = T)
data_6 <- read.csv("data/Fig5g-otu-heatmap-2-ah-MH-AA.csv", header = T)
data_7 <- read.csv("data/Fig5g-otu-heatmap-2-ah-WYJ-RA.csv", header = T)
data_8 <- read.csv("data/Fig5g-otu-heatmap-2-ah-WYJ-AA.csv", header = T)

all_1 <- gather(data_1, key = group , value = logFC, HnMH63RA)
all_2 <- gather(data_2, key = group , value = logFC, HnMH63AA)
all_3 <- gather(data_3, key = group , value = logFC, HnWYJRA)
all_4 <- gather(data_4, key = group , value = logFC, HnWYJAA)
all_5 <- gather(data_5, key = group , value = logFC, AhMH63RA)
all_6 <- gather(data_6, key = group , value = logFC, AhMH63AA)
all_7 <- gather(data_7, key = group , value = logFC, AhWYJRA)
all_8 <- gather(data_8, key = group , value = logFC, AhWYJAA)



p1 <- ggplot(all_1, aes(group, HnMH63)) + 
  geom_tile(aes(fill = logFC)) + 
  scale_fill_gradient2(low = "#6D9EC1",high = "sandybrown", midpoint = 0)+
  theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust = 0))
p1

p2 <- ggplot(all_2, aes(group, HnMH63)) + 
  geom_tile(aes(fill = logFC)) + 
  scale_fill_gradient2(low = "#6D9EC1",high = "sandybrown", midpoint = 0)+
  theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust = 0))
p2


p3 <- ggplot(all_3, aes(group, HnWYJ)) + 
  geom_tile(aes(fill = logFC)) + 
  scale_fill_gradient2(low = "#6D9EC1",high = "sandybrown", midpoint = 0)+
  theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust = 0))
p3

p4 <- ggplot(all_4, aes(group, HnWYJ)) + 
  geom_tile(aes(fill = logFC)) + 
  scale_fill_gradient2(low = "#6D9EC1",high = "sandybrown", midpoint = 0)+
  theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust = 0))
p4


p5 <- ggplot(all_5, aes(group, AhMH63)) + 
  geom_tile(aes(fill = logFC)) + 
  scale_fill_gradient2(low = "#6D9EC1",high = "sandybrown", midpoint = 0)+
  theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust = 0))
p5

p6 <- ggplot(all_6, aes(group, AhMH63)) + 
  geom_tile(aes(fill = logFC)) + 
  scale_fill_gradient2(low = "#6D9EC1",high = "sandybrown", midpoint = 0)+
  theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust = 0))
p6


p7 <- ggplot(all_7, aes(group, AhWYJ)) + 
  geom_tile(aes(fill = logFC)) + 
  scale_fill_gradient2(low = "#6D9EC1",high = "sandybrown", midpoint = 0)+
  theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust = 0))
p7

p8 <- ggplot(all_8, aes(group, AhWYJ)) + 
  geom_tile(aes(fill = logFC)) + 
  scale_fill_gradient2(low = "#6D9EC1",high = "sandybrown", midpoint = 0)+
  theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust = 0))
p8


ggsave(paste("result/Fig5g-Heatmap-HN-1.pdf", sep=""), p1, width = 2, height = 40)
ggsave(paste("result/Fig5g-Heatmap-HN-2.pdf", sep=""), p2, width = 2, height = 40)
ggsave(paste("result/Fig5g-Heatmap-HN-3.pdf", sep=""), p3, width = 2, height = 40)
ggsave(paste("result/Fig5g-Heatmap-HN-4.pdf", sep=""), p4, width = 2, height = 40)

ggsave(paste("result/Fig5g-Heatmap-AH-1.pdf", sep=""), p5, width = 2, height = 40)
ggsave(paste("result/Fig5g-Heatmap-AH-2.pdf", sep=""), p6, width = 2, height = 40)
ggsave(paste("result/Fig5g-Heatmap-AH-3.pdf", sep=""), p7, width = 2, height = 40)
ggsave(paste("result/Fig5g-Heatmap-AH-4.pdf", sep=""), p8, width = 2, height = 40)



