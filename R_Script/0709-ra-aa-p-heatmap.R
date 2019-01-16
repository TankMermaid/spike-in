setwd("/mnt/bai/qinyuan/xiaoxuan/bac_New/")
print(paste("Your working directory is in",getwd()))

library(ggplot2)
library(plyr)
library(scales)
library(reshape2)

## All

data <- read.csv("table/20180724-heatmap-p.csv", header = T)

## HN 12个门

data_hn <- select(data, Bacteria, starts_with("HN"))

all <- gather(data_hn, key = group , value = logFC, `HN_MH_RA`:`HN_WYJ_AA`)

p_hn <- ggplot(all, aes(group, Bacteria)) + 
  geom_tile(aes(fill = logFC),colour = "grey50") + 
  scale_fill_gradient2(low = "#6D9EC1",high = "sandybrown", midpoint = 0)+
  theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust = 0))+
  scale_y_discrete(limits=c("Gammaproteobacteria","Deltaproteobacteria","Betaproteobacteria",
                            "Alphaproteobacteria","Verrucomicrobia","Spirochaetes","Nitrospirae",
                            "Ignavibacteriae","Firmicutes","Chloroflexi","Bacteroidetes","Actinobacteria",
                            "Acidobacteria"))+
  scale_x_discrete(limits=c("HN_MH_RA","HN_MH_AA","HN_WYJ_RA","HN_WYJ_AA"))
p_hn

ggsave(paste("figure/20180725-Heatmap-HN.pdf", sep=""), p_hn, width = 3.5, height = 4)

### 不筛选 18个门

p_hn_18 <- ggplot(all, aes(group, Bacteria)) + 
  geom_tile(aes(fill = logFC),colour = "grey50") + 
  scale_fill_gradient2(low = "#6D9EC1",high = "sandybrown", midpoint = 0)+
  theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust = 0))+
  scale_y_discrete(limits=c("Gammaproteobacteria","Deltaproteobacteria","Betaproteobacteria",
                            "Alphaproteobacteria","Verrucomicrobia","Tenericutes","Spirochaetes","Poribacteria","Planctomycetes","Nitrospirae","Latescibacteria",
                            "Ignavibacteriae","Gemmatimonadetes","Firmicutes","Chloroflexi","Bacteroidetes","Actinobacteria",
                            "Acidobacteria"))+
  scale_x_discrete(limits=c("HN_MH_RA","HN_MH_AA","HN_WYJ_RA","HN_WYJ_AA"))
p_hn_18

ggsave(paste("figure/20180725-Heatmap-HN-18.pdf", sep=""), p_hn_18, width = 3.5, height = 4)


## AH 12个门

data_ah <- select(data, Bacteria, starts_with("AH"))

all2 <- gather(data_ah, key = group , value = logFC, `AH_MH_RA`:`AH_WYJ_AA`)

p_ah <- ggplot(all2, aes(group, Bacteria)) + 
  geom_tile(aes(fill = logFC),colour = "grey50") + 
  scale_fill_gradient2(low = "#6D9EC1",high = "sandybrown", midpoint = 0)+
  theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust = 0))+
  scale_y_discrete(limits=c("Gammaproteobacteria","Deltaproteobacteria","Betaproteobacteria",
                            "Alphaproteobacteria","Verrucomicrobia","Spirochaetes","Nitrospirae",
                            "Firmicutes","Chloroflexi","Bacteroidetes","Actinobacteria",
                            "Acidobacteria"))+
  scale_x_discrete(limits=c("AH_MH_RA","AH_MH_AA","AH_WYJ_RA","AH_WYJ_AA"))
p_ah

ggsave(paste("figure/20180725-Heatmap-AH.pdf", sep=""), p_ah, width = 3.5, height = 4)

### 不筛选 18个门

p_ah_18 <- ggplot(all2, aes(group, Bacteria)) + 
  geom_tile(aes(fill = logFC),colour = "grey50") + 
  scale_fill_gradient2(low = "#6D9EC1",high = "sandybrown", midpoint = 0)+
  theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust = 0))+
  scale_y_discrete(limits=c("Gammaproteobacteria","Deltaproteobacteria","Betaproteobacteria",
                            "Alphaproteobacteria","Verrucomicrobia","Tenericutes","Spirochaetes","Poribacteria","Planctomycetes","Nitrospirae","Latescibacteria",
                            "Ignavibacteriae","Gemmatimonadetes","Firmicutes","Chloroflexi","Bacteroidetes","Actinobacteria",
                            "Acidobacteria"))+
  scale_x_discrete(limits=c("AH_MH_RA","AH_MH_AA","AH_WYJ_RA","AH_WYJ_AA"))
p_ah_18

ggsave(paste("figure/20180725-Heatmap-AH-18.pdf", sep=""), p_ah_18, width = 3.5, height = 4)




## 不算logFC ra


data2 <- read.csv("table/20180720-heatmap-ra.csv", header = T)
all2 <- gather(data2, key = group , value = logFC, `AH_MH_D`:`HN_WYJ_W`)

p2 <- ggplot(all2, aes(group, Bacteria)) + 
  geom_tile(aes(fill = logFC),colour = "grey50") + 
  scale_fill_gradient2(low = "#6D9EC1",high = "sandybrown", midpoint = 0)+
  theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust = 0))+
  scale_x_discrete(limits=c("AH_MH_D","AH_MH_W","HN_MH_D","HN_MH_W","AH_WYJ_D","AH_WYJ_W","HN_WYJ_D",
                            "HN_WYJ_W"))
p2


ggsave(paste("figure/20180717-Heatmap-All4.pdf", sep=""), p2, width = 7, height = 5)

## 不算logFC AA

data3 <- read.csv("table/20180720-heatmap-aa.csv", header = T)
all3 <- gather(data3, key = group , value = logFC, `AH_MH_D`:`HN_WYJ_W`)

p3 <- ggplot(all3, aes(group, Bacteria)) + 
  geom_tile(aes(fill = logFC),colour = "grey50") + 
  scale_fill_gradient2(low = "#6D9EC1",high = "sandybrown", midpoint = 0)+
  theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust = 0))+
  scale_x_discrete(limits=c("AH_MH_D","AH_MH_W","HN_MH_D","HN_MH_W","AH_WYJ_D","AH_WYJ_W","HN_WYJ_D",
                            "HN_WYJ_W"))
p3



ggsave(paste("figure/20180717-Heatmap-All5.pdf", sep=""), p3, width = 7, height = 5)

