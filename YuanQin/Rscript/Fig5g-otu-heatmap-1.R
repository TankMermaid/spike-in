setwd("/mnt/bai/qinyuan/xiaoxuan/QA/")
print(paste("Your working directory is in",getwd()))

library(ggplot2)
library(plyr)
library(scales)
library(reshape2)
library(pheatmap)


data_hn <- read.csv("data/Fig5g-otu-heatmap-1-hn.csv", header = T)
data_ah <- read.csv("data/Fig5g-otu-heatmap-1-ah.csv", header = T)

## HN

all_hn <- gather(data_hn, key = group , value = logFC, `HnMH63_RA`:`HnWYJ_AA`)

p_hn <- ggplot(all_hn, aes(group, OTUS)) + 
  geom_tile(aes(fill = logFC)) + 
  scale_fill_gradient2(low = "#6D9EC1",high = "sandybrown", midpoint = 0)+
  theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust = 0))+
  scale_x_discrete(limits=c("HnMH63_RA","HnMH63_AA","HnWYJ_RA","HnWYJ_AA"))
p_hn

## AH

all_ah <- gather(data_ah, key = group , value = logFC, `AhMH63_RA`:`AhWYJ_AA`)

p_ah <- ggplot(all_ah, aes(group, OTUS)) + 
  geom_tile(aes(fill = logFC)) + 
  scale_fill_gradient2(low = "#6D9EC1",high = "sandybrown", midpoint = 0)+
  theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust = 0))+
  scale_x_discrete(limits=c("AhMH63_RA","AhMH63_AA","AhWYJ_RA","AhWYJ_AA"))
p_ah

ggsave(paste("result/Fig5g-Heatmap-HNALL.pdf", sep=""), p_hn, width = 3, height = 40)
ggsave(paste("result/Fig5g-Heatmap-AHALL.pdf", sep=""), p_ah, width = 3, height = 40)


#############################################################################
## HN

data_hn <- select(data, OTUS, starts_with("Hn"))
data_hn1 <- data.frame(data_hn[1:140,])
data_hn2 <- data.frame(data_hn[141:280,])
data_hn3 <- data.frame(data_hn[281:420,])
data_hn4 <- data.frame(data_hn[421:560,])

all1 <- gather(data_hn1, key = group , value = logFC, `HnMH63_RA`:`HnWYJ_AA`)
all2 <- gather(data_hn2, key = group , value = logFC, `HnMH63_RA`:`HnWYJ_AA`)
all3 <- gather(data_hn3, key = group , value = logFC, `HnMH63_RA`:`HnWYJ_AA`)
all4 <- gather(data_hn4, key = group , value = logFC, `HnMH63_RA`:`HnWYJ_AA`)

p_hn1 <- ggplot(all1, aes(group, OTUS)) + 
  geom_tile(aes(fill = logFC)) + 
  scale_fill_gradientn(colours=rainbow(7))+
  theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust = 0))+
  scale_x_discrete(limits=c("HnMH63_RA","HnMH63_AA","HnWYJ_RA","HnWYJ_AA"))
p_hn1

p_hn2 <- ggplot(all2, aes(group, OTUS)) + 
  geom_tile(aes(fill = logFC)) + 
  scale_fill_gradientn(colours=rainbow(7))+
  theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust = 0))+
  scale_x_discrete(limits=c("HnMH63_RA","HnMH63_AA","HnWYJ_RA","HnWYJ_AA"))
p_hn2

p_hn3 <- ggplot(all3, aes(group, OTUS)) + 
  geom_tile(aes(fill = logFC)) + 
  scale_fill_gradientn(colours=rainbow(7))+
  theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust = 0))+
  scale_x_discrete(limits=c("HnMH63_RA","HnMH63_AA","HnWYJ_RA","HnWYJ_AA"))
p_hn3


p_hn4 <- ggplot(all4, aes(group, OTUS)) + 
  geom_tile(aes(fill = logFC)) + 
  scale_fill_gradientn(colours=rainbow(7))+
  theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust = 0))+
  scale_x_discrete(limits=c("HnMH63_RA","HnMH63_AA","HnWYJ_RA","HnWYJ_AA"))
p_hn4

ggsave(paste("result/Fig5g-Heatmap-HN1.pdf", sep=""), p_hn1, width = 3, height = 9)
ggsave(paste("result/Fig5g-Heatmap-HN2.pdf", sep=""), p_hn2, width = 3, height = 9)
ggsave(paste("result/Fig5g-Heatmap-HN3.pdf", sep=""), p_hn3, width = 3, height = 9)
ggsave(paste("result/Fig5g-Heatmap-HN4.pdf", sep=""), p_hn4, width = 3, height = 9)



## AH

data_ah <- select(data, OTUS, starts_with("Ah"))
data_ah1 <- data.frame(data_ah[1:140,])
data_ah2 <- data.frame(data_ah[141:280,])
data_ah3 <- data.frame(data_ah[281:420,])
data_ah4 <- data.frame(data_ah[421:560,])

all1 <- gather(data_ah1, key = group , value = logFC, `AhMH63_RA`:`AhWYJ_AA`)
all2 <- gather(data_ah2, key = group , value = logFC, `AhMH63_RA`:`AhWYJ_AA`)
all3 <- gather(data_ah3, key = group , value = logFC, `AhMH63_RA`:`AhWYJ_AA`)
all4 <- gather(data_ah4, key = group , value = logFC, `AhMH63_RA`:`AhWYJ_AA`)

p_ah1 <- ggplot(all1, aes(group, OTUS)) + 
  geom_tile(aes(fill = logFC)) + 
  scale_fill_gradientn(colours=rainbow(7))+
  theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust = 0))+
  scale_x_discrete(limits=c("AhMH63_RA","AhMH63_AA","AhWYJ_RA","AhWYJ_AA"))
p_ah1

p_ah2 <- ggplot(all2, aes(group, OTUS)) + 
  geom_tile(aes(fill = logFC)) + 
  scale_fill_gradientn(colours=rainbow(7))+
  theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust = 0))+
  scale_x_discrete(limits=c("AhMH63_RA","AhMH63_AA","AhWYJ_RA","AhWYJ_AA"))
p_ah2

p_ah3 <- ggplot(all3, aes(group, OTUS)) + 
  geom_tile(aes(fill = logFC)) + 
  scale_fill_gradientn(colours=rainbow(7))+
  theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust = 0))+
  scale_x_discrete(limits=c("AhMH63_RA","AhMH63_AA","AhWYJ_RA","AhWYJ_AA"))
p_ah3


p_ah4 <- ggplot(all4, aes(group, OTUS)) + 
  geom_tile(aes(fill = logFC)) + 
  scale_fill_gradientn(colours=rainbow(7))+
  theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust = 0))+
  scale_x_discrete(limits=c("AhMH63_RA","AhMH63_AA","AhWYJ_RA","AhWYJ_AA"))
p_ah4

ggsave(paste("result/Fig5g-Heatmap-ah1.pdf", sep=""), p_ah1, width = 3, height = 9)
ggsave(paste("result/Fig5g-Heatmap-ah2.pdf", sep=""), p_ah2, width = 3, height = 9)
ggsave(paste("result/Fig5g-Heatmap-ah3.pdf", sep=""), p_ah3, width = 3, height = 9)
ggsave(paste("result/Fig5g-Heatmap-ah4.pdf", sep=""), p_ah4, width = 3, height = 9)
